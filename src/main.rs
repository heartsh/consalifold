extern crate consalifold;
extern crate num_cpus;
extern crate crossbeam;

use consalifold::*;
use std::env;
use std::path::Path;
use std::io::{BufRead, BufWriter};
use std::fs::File;
use std::fs::create_dir;
use crossbeam::scope;
use std::process::{Command, Output};
use std::fs::remove_file;

type MeaCssStr = MeaSsStr;

const DEFAULT_MIX_WEIGHT: Prob = 0.5;
const MIN_POW_OF_2: i32 = -4;
const MAX_POW_OF_2: i32 = 10;
const DEFAULT_GAMMA: Prob = NEG_INFINITY;
enum ScoringModel {
  Turner,
  Contra,
  Posterior,
}
const DEFAULT_SCORING_MODEL: &str = "turner";
const README_CONTENTS_2: &str = "# gamma=x.sth\nThis file type contains a predicted consensus secondary structure in Stockholm format, and this predicted consensus structure is under the prediction accuracy control parameter \"x.\"\n\n";

fn main() {
  let args = env::args().collect::<Vec<Arg>>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "A path to an input CLUSTAL/FASTA/STOCKHOLM file containing the sequence alignment of RNA sequences", "STR");
  opts.reqopt("o", "output_dir_path", "The path to an output directory", "STR");
  opts.optopt("", "min_base_pair_prob", &format!("A minimum base-pairing-probability (Uses {} (Turner)/{}(CONTRAfold) by default)", DEFAULT_MIN_BPP, DEFAULT_MIN_BPP_CONTRA), "FLOAT");
  opts.optopt("", "offset_4_max_gap_num", &format!("An offset for maximum numbers of gaps (Uses {} by default)", DEFAULT_OFFSET_4_MAX_GAP_NUM), "UINT");
  opts.optopt("g", "gamma", "A specific gamma parameter rather than a range of gamma parameters", "FLOAT");
  opts.optopt("", "mix_weight", &format!("A mixture weight (Uses {} by default)", DEFAULT_MIX_WEIGHT), "FLOAT");
  opts.optopt("m", "scoring_model", &format!("Choose a structural alignment scoring model from turner, contra, posterior (Uses {} by default)", DEFAULT_SCORING_MODEL), "STR");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of all the threads of this computer by default)", "UINT");
  opts.optflag("s", "produces_struct_profs", &format!("Also compute RNA structural context profiles"));
  opts.optflag("p", "outputs_probs", &format!("Output probabilities"));
  opts.optflag("h", "help", "Print a help menu");
  let matches = match opts.parse(&args[1 ..]) {
    Ok(opt) => {opt}
    Err(failure) => {print_program_usage(&program_name, &opts); panic!(failure.to_string())}
  };
  if matches.opt_present("h") {
    print_program_usage(&program_name, &opts);
    return;
  }
  let input_file_path = matches.opt_str("i").unwrap();
  let input_file_path = Path::new(&input_file_path);
  let output_dir_path = matches.opt_str("o").unwrap();
  let output_dir_path = Path::new(&output_dir_path);
  let scoring_model = if matches.opt_present("m") {
    let scoring_model_str = matches.opt_str("m").unwrap();
    if scoring_model_str == "turner" {
      ScoringModel::Turner
    } else if scoring_model_str == "contra" {
      ScoringModel::Contra
    } else if scoring_model_str == "posterior" {
      ScoringModel::Posterior
    } else {
      assert!(false);
      ScoringModel::Turner
    }
  } else {
    ScoringModel::Turner
  };
  let uses_contra_model = matches!(scoring_model, ScoringModel::Contra);
  let min_bpp = if matches.opt_present("min_base_pair_prob") {
    matches.opt_str("min_base_pair_prob").unwrap().parse().unwrap()
  } else {
    if uses_contra_model {DEFAULT_MIN_BPP_CONTRA} else {DEFAULT_MIN_BPP}
  };
  let offset_4_max_gap_num = if matches.opt_present("offset_4_max_gap_num") {
    matches.opt_str("offset_4_max_gap_num").unwrap().parse().unwrap()
  } else {
    DEFAULT_OFFSET_4_MAX_GAP_NUM
  };
  let gamma = if matches.opt_present("gamma") {
    matches.opt_str("gamma").unwrap().parse().unwrap()
  } else {
    DEFAULT_GAMMA
  };
  let mix_weight = if matches.opt_present("mix_weight") {
    matches.opt_str("mix_weight").unwrap().parse().unwrap()
  } else {
    DEFAULT_MIX_WEIGHT
  };
  let produces_struct_profs = matches.opt_present("s");
  let outputs_probs = matches.opt_present("p") || produces_struct_profs;
  let num_of_threads = if matches.opt_present("t") {
    matches.opt_str("t").unwrap().parse().unwrap()
  } else {
    num_cpus::get() as NumOfThreads
  };
  let extension = input_file_path.extension().unwrap().to_str().unwrap();
  let is_stockholm = match extension {
    "sto" | "stk" | "sth" => true,
    _ => false,
  };
  let is_fasta = match extension {
    "fasta" | "fna" | "ffn" | "faa" | "frn" | "fa" => true,
    _ => false,
  };
  let (cols, seq_ids) = if is_stockholm {
    read_sa_from_stockholm_file(input_file_path)
  } else if is_fasta {
    read_sa_from_fasta_file(input_file_path)
  } else {
    read_sa_from_clustal_file(input_file_path)
  };
  let sa_len = cols.len();
  let mut thread_pool = Pool::new(num_of_threads);
  if sa_len + 2 <= u8::MAX as usize {
    multi_threaded_consalifold::<u8>(&mut thread_pool, &cols, &seq_ids, offset_4_max_gap_num, min_bpp, produces_struct_profs, output_dir_path, gamma, outputs_probs, mix_weight, input_file_path, scoring_model);
  } else {
    multi_threaded_consalifold::<u16>(&mut thread_pool, &cols, &seq_ids, offset_4_max_gap_num, min_bpp, produces_struct_profs, output_dir_path, gamma, outputs_probs, mix_weight, input_file_path, scoring_model);
  }
}

fn multi_threaded_consalifold<T>(thread_pool: &mut Pool, cols: &Cols, seq_ids: &SeqIds, offset_4_max_gap_num: usize, min_bpp: Prob, produces_struct_profs: bool, output_dir_path: &Path, gamma: Prob, outputs_probs: bool, mix_weight: Prob, input_file_path: &Path, scoring_model: ScoringModel)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Display + Sync + Send,
{
  let uses_contra_model = matches!(scoring_model, ScoringModel::Contra);
  let is_posterior_model = matches!(scoring_model, ScoringModel::Posterior);
  let mut sa = SeqAlign::<T>::new();
  sa.cols = cols.clone();
  let num_of_rnas = sa.cols[0].len();
  let mut seq_lens = vec![0 as usize; num_of_rnas];
  let sa_len = sa.cols.len();
  sa.pos_map_sets = vec![vec![T::zero(); num_of_rnas]; sa_len];
  let mut fasta_records = vec![FastaRecord::origin(); num_of_rnas];
  for i in 0 .. sa_len {
    for j in 0 .. num_of_rnas {
      let base = sa.cols[i][j];
      if base != PSEUDO_BASE {
        fasta_records[j].seq.push(base);
        seq_lens[j] += 1;
        sa.pos_map_sets[i][j] = T::from_usize(seq_lens[j]).unwrap();
      }
    }
  }
  for i in 0 .. num_of_rnas {
    fasta_records[i].seq.insert(0, PSEUDO_BASE);
    fasta_records[i].seq.push(PSEUDO_BASE);
    fasta_records[i].fasta_id = seq_ids[i].clone();
  }
  let ref ref_2_sa = sa;
  let ref ref_2_fasta_records = fasta_records;
  let mut mix_bpp_mat = ProbMat::new();
  let ref mut ref_2_mix_bpp_mat = mix_bpp_mat;
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  scope(|scope| {
    let handler = scope.spawn(|_| {
      let output_file_prefix = input_file_path.file_stem().unwrap().to_str().unwrap();
      let arg = format!("--id-prefix={}", output_file_prefix);
      let args = vec!["-p", input_file_path.to_str().unwrap(), &arg, "--noPS", "--noDP"];
      let _ = run_command("RNAalifold", &args, "Failed to run RNAalifold");
      let mut rnaalifold_bpp_mat = SparseProbMat::<T>::default();
      let cwd = env::current_dir().unwrap();
      let output_file_path = cwd.join(String::from(output_file_prefix) + "_0001_ali.out");
      let output_file = BufReader::new(File::open(output_file_path.clone()).unwrap());
      for (k, line) in output_file.lines().enumerate() {
        if k == 0 {continue;}
        let line = line.unwrap();
        if !line.starts_with(" ") {continue;}
        let substrings: Vec<&str> = line.split_whitespace().collect();
        let i = T::from_usize(substrings[0].parse().unwrap()).unwrap() - T::one();
        let j = T::from_usize(substrings[1].parse().unwrap()).unwrap() - T::one();
        let mut bpp = String::from(substrings[3]);
        bpp.pop();
        let bpp = 0.01 * bpp.parse::<Prob>().unwrap();
        if bpp == 0. {continue;}
        rnaalifold_bpp_mat.insert((i, j), bpp);
      }
      let _ = remove_file(output_file_path);
      rnaalifold_bpp_mat
    });
    let prob_mat_sets = if !is_posterior_model {consprob::<T>(thread_pool, ref_2_fasta_records, min_bpp, T::from_usize(offset_4_max_gap_num).unwrap(), produces_struct_profs, uses_contra_model)} else {locarnap_plus_pct(thread_pool, ref_2_fasta_records, output_dir_path)};
    if outputs_probs {
      write_prob_mat_sets::<T>(output_dir_path, &prob_mat_sets, produces_struct_profs);
    }
    let rnaalifold_bpp_mat = handler.join().unwrap();
    *ref_2_mix_bpp_mat = get_mix_bpp_mat(&prob_mat_sets, &rnaalifold_bpp_mat, ref_2_sa, mix_weight);
  }).unwrap();
  if gamma != NEG_INFINITY {
    let output_file_path = output_dir_path.join(&format!("gamma={}.sth", gamma));
    compute_and_write_mea_css(&mix_bpp_mat, &sa, gamma, &output_file_path, &fasta_records);
  } else {
    thread_pool.scoped(|scope| {
      for pow_of_2 in MIN_POW_OF_2 .. MAX_POW_OF_2 + 1 {
        let gamma = (2. as Prob).powi(pow_of_2);
        let ref ref_2_mix_bpp_mat = mix_bpp_mat;
        let ref ref_2_sa = sa;
        let ref ref_2_fasta_records = fasta_records;
        let output_file_path = output_dir_path.join(&format!("gamma={}.sth", gamma));
        scope.execute(move || {
          compute_and_write_mea_css::<T>(ref_2_mix_bpp_mat, ref_2_sa, gamma, &output_file_path, ref_2_fasta_records);
        });
      }
    });
  }
  let mut readme_contents = String::from(README_CONTENTS_2);
  readme_contents.push_str(README_CONTENTS);
  write_readme(output_dir_path, &readme_contents);
}

fn locarnap_plus_pct<T>(thread_pool: &mut Pool, fasta_records: &FastaRecords, output_dir_path: &Path) -> ProbMatSets<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Sync + Send,
{
  let num_of_fasta_records = fasta_records.len();
  for i in 0 .. num_of_fasta_records {
    let output_file_path = output_dir_path.join(&format!("locarnap_seq_{}.fa", i));
    let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).unwrap());
    let ref seq = fasta_records[i].seq;
    let seq_len = seq.len();
    let buf_4_writer_2_output_file = format!(">{}\n{}", i, revert(&seq[1 .. seq_len - 1]));
    let _ = writer_2_output_file.write_all(buf_4_writer_2_output_file.as_bytes());
  }
  let mut prob_mats_with_rna_id_pairs = StaProbMatsWithRnaIdPairs::<T>::default();
  for rna_id_1 in 0 .. num_of_fasta_records {
    for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
      let rna_id_pair = (rna_id_1, rna_id_2);
      prob_mats_with_rna_id_pairs.insert(rna_id_pair, StaProbMats::<T>::origin());
    }
  }
  thread_pool.scoped(|scope| {
    for (rna_id_pair, prob_mats) in prob_mats_with_rna_id_pairs.iter_mut() {
      scope.execute(move || {
        prob_mats.bpp_mat_pair = exec_locarnap(rna_id_pair, output_dir_path);
      });
    }
  });
  let mut prob_mat_sets = vec![PctStaProbMats::<T>::origin(); num_of_fasta_records];
  thread_pool.scoped(|scope| {
    for (rna_id, prob_mats) in prob_mat_sets.iter_mut().enumerate() {
      let ref ref_2_prob_mats_with_rna_id_pairs = prob_mats_with_rna_id_pairs;
      let seq_len = fasta_records[rna_id].seq.len();
      let output_file_path = output_dir_path.join(&format!("locarnap_seq_{}.fa", rna_id));
      scope.execute(move || {
        *prob_mats = pct_of_bpp_mats_locarnap::<T>(ref_2_prob_mats_with_rna_id_pairs, rna_id, seq_len, num_of_fasta_records);
        let _ = remove_file(output_file_path);
      });
    }
  });
  prob_mat_sets
}

fn exec_locarnap<T>(rna_id_pair: &RnaIdPair, output_dir_path: &Path) -> SparseProbMatPair<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Sync + Send,
{
  let mut bpp_mat_pair = (
    SparseProbMat::<T>::default(),
    SparseProbMat::<T>::default(),
    );
  let (seq_file_path_1, seq_file_path_2) = (
    output_dir_path.join(&format!("locarnap_seq_{}.fa", rna_id_pair.0)),
    output_dir_path.join(&format!("locarnap_seq_{}.fa", rna_id_pair.1)),
    );
  let output_file_path = output_dir_path.join(format!("locarnap_seq_{}_seq_{}.dat", rna_id_pair.0, rna_id_pair.1));
  let arg = format!("--write-arcmatch-probs={}", output_file_path.to_str().unwrap());
  let args = vec![seq_file_path_1.to_str().unwrap(), seq_file_path_2.to_str().unwrap(), &arg];
  let _ = run_command("locarna_p", &args, "Failed to run LocARNA-P");
  let output_file = BufReader::new(File::open(output_file_path.clone()).unwrap());
  for line in output_file.lines() {
    let strs: Vec<String> = line.unwrap().trim().split_whitespace().map(|x| String::from(x)).collect();
    let (i, j, k, l, bpap) = (
      T::from_usize(strs[0].parse().unwrap()).unwrap(),
      T::from_usize(strs[1].parse().unwrap()).unwrap(),
      T::from_usize(strs[2].parse().unwrap()).unwrap(),
      T::from_usize(strs[3].parse().unwrap()).unwrap(),
      strs[4].parse().unwrap(),
      );
    match bpp_mat_pair.0.get_mut(&(i, j)) {
      Some(bpp) => {
        *bpp += bpap;
      }, None => {
        bpp_mat_pair.0.insert((i, j), bpap);
      },
    }
    match bpp_mat_pair.1.get_mut(&(k, l)) {
      Some(bpp) => {
        *bpp += bpap;
      }, None => {
        bpp_mat_pair.1.insert((k, l), bpap);
      },
    }
  }
  let _ = remove_file(output_file_path);
  bpp_mat_pair
}

pub fn revert<'a>(seq: &'a [usize]) -> String {
  let mut new_seq = Vec::<u8>::new();
  for &c in seq {
    let new_base = match c {
      A => BIG_A as u8,
      C => BIG_C as u8,
      G => BIG_G as u8,
      U => BIG_U as u8,
      _ => {assert!(false); U as u8},
    };
    new_seq.push(new_base);
  }
  String::from_utf8(new_seq).unwrap()
}

fn pct_of_bpp_mats_locarnap<T>(prob_mats_with_rna_id_pairs: &StaProbMatsWithRnaIdPairs<T>, rna_id: RnaId, upp_mat_len: usize, num_of_rnas: usize, ) -> PctStaProbMats<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let weight = 1. / (num_of_rnas - 1) as Prob;
  let mut pct_prob_mats = PctStaProbMats::new(upp_mat_len);
  for rna_id_2 in 0 .. num_of_rnas {
    if rna_id == rna_id_2 {continue;}
    let rna_id_pair = if rna_id < rna_id_2 {(rna_id, rna_id_2)} else {(rna_id_2, rna_id)};
    let ref ref_2_prob_mats = prob_mats_with_rna_id_pairs[&rna_id_pair];
    let ref_2_bpp_mat = if rna_id < rna_id_2 {
      &ref_2_prob_mats.bpp_mat_pair.0
    } else {
      &ref_2_prob_mats.bpp_mat_pair.1
    };
    for (pos_pair, &bpp) in ref_2_bpp_mat.iter() {
      let weighted_bpp = weight * bpp;
      match pct_prob_mats.bpp_mat.get_mut(pos_pair) {
        Some(bpp) => {
          *bpp += weighted_bpp;
        },
        None => {
          pct_prob_mats.bpp_mat.insert(*pos_pair, weighted_bpp);
        },
      }
    }
  }
  pct_prob_mats
}

fn get_mix_bpp_mat<T>(prob_mat_sets: &ProbMatSets<T>, rnaalifold_bpp_mat: &SparseProbMat<T>, sa: &SeqAlign<T>, mix_weight: Prob) -> ProbMat
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Display + Sync + Send,
{
  let sa_len = sa.cols.len();
  let num_of_rnas = sa.cols[0].len();
  let mut mix_bpp_mat = vec![vec![0.; sa_len]; sa_len];
  for i in 0 .. sa_len {
    for j in i + 1 .. sa_len {
      let mut mean_bpp = 0.;
      let mut effective_num_of_rnas = 0;
      for k in 0 .. num_of_rnas {
        if sa.cols[i][k] == PSEUDO_BASE || sa.cols[j][k] == PSEUDO_BASE {continue;}
        let ref bpp_mat = prob_mat_sets[k].bpp_mat;
        let pos_pair = (sa.pos_map_sets[i][k], sa.pos_map_sets[j][k]);
        match bpp_mat.get(&pos_pair) {
          Some(&bpp) => {
            mean_bpp += bpp;
            effective_num_of_rnas += 1;
          }, None => {},
        }
      }
      mix_bpp_mat[i][j] = if effective_num_of_rnas > 0 {mix_weight * mean_bpp / effective_num_of_rnas as Prob} else {0.};
      let pos_pair = (T::from_usize(i).unwrap(), T::from_usize(j).unwrap());
      match rnaalifold_bpp_mat.get(&pos_pair) {
        Some(&rnaalifold_bpp) => {
          mix_bpp_mat[i][j] += (1. - mix_weight) * rnaalifold_bpp;
        }, None => {},
      }
    }
  }
  mix_bpp_mat
}

fn compute_and_write_mea_css<T>(mix_bpp_mat: &ProbMat, sa: &SeqAlign<T>, gamma: Prob, output_file_path: &Path, fasta_records: &FastaRecords)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Display + Sync + Send,
{
  let mea_css = consalifold::<T>(mix_bpp_mat, gamma, sa);
  let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).unwrap());
  let mut buf_4_writer_2_output_file = format!("# STOCKHOLM 1.0\n");
  let sa_len = sa.cols.len();
  let max_seq_id_len = fasta_records.iter().map(|fasta_record| {fasta_record.fasta_id.len()}).max().unwrap();
  let num_of_rnas = sa.cols[0].len();
  for rna_id in 0 .. num_of_rnas {
    let ref seq_id = fasta_records[rna_id].fasta_id;
    buf_4_writer_2_output_file.push_str(seq_id);
    let mut stockholm_row = vec![' ' as Char; max_seq_id_len - seq_id.len() + 2];
    let mut sa_row = (0 .. sa_len).map(|x| {revert_char(sa.cols[x][rna_id])}).collect::<Vec<Char>>();
    stockholm_row.append(&mut sa_row);
    let stockholm_row = unsafe {from_utf8_unchecked(&stockholm_row)};
    buf_4_writer_2_output_file.push_str(&stockholm_row);
    buf_4_writer_2_output_file.push_str("\n");
  }
  let descriptor = "#=GC SS_cons";
  let descriptor_len = descriptor.len();
  buf_4_writer_2_output_file.push_str(descriptor);
  let mut stockholm_row = vec![' ' as Char; max_seq_id_len - descriptor_len + 2];
  let mut mea_css_str = get_mea_css_str(&mea_css, sa_len);
  stockholm_row.append(&mut mea_css_str);
  let stockholm_row = unsafe {from_utf8_unchecked(&stockholm_row)};
  buf_4_writer_2_output_file.push_str(&stockholm_row);
  buf_4_writer_2_output_file.push_str("\n//");
  let _ = writer_2_output_file.write_all(buf_4_writer_2_output_file.as_bytes());
}

fn get_mea_css_str<T>(mea_css: &MeaCss<T>, sa_len: usize) -> MeaCssStr
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Display + Sync + Send,
{
  let mut mea_css_str = vec![UNPAIRING_BASE; sa_len];
  for &(i, j) in &mea_css.bpa_pos_pairs {
    mea_css_str[i.to_usize().unwrap()] = BASE_PAIRING_LEFT_BASE;
    mea_css_str[j.to_usize().unwrap()] = BASE_PAIRING_RIGHT_BASE;
  }
  mea_css_str
}

fn run_command(command: &str, args: &[&str], expect: &str) -> Output {
  Command::new(command).args(args).output().expect(expect)
}
