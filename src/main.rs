extern crate consalifold;
extern crate num_cpus;

use consalifold::*;
use std::env;
use std::path::Path;
use std::io::{BufReader, BufWriter};
use std::fs::File;
use std::fs::create_dir;

type MeaCssStr = MeaSsStr;

const DEFAULT_MIX_WEIGHT: Prob = 0.5;
const DEFAULT_MIN_POW_OF_2: i32 = -4;
const DEFAULT_MAX_POW_OF_2: i32 = 10;
const GAMMA_4_BENCH: Prob = 1.;

fn main() {
  let args = env::args().collect::<Vec<Arg>>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "A path to an input CLUSTAL file containing the sequence alignment of RNA sequences", "STR");
  opts.reqopt("o", "output_dir_path", "The path to an output directory", "STR");
  opts.optopt("", "min_base_pair_prob", &format!("A minimum base-pairing-probability (Uses {} by default)", DEFAULT_MIN_BPP), "FLOAT");
  opts.optopt("", "offset_4_max_gap_num", &format!("An offset for maximum numbers of gaps (Uses {} by default)", DEFAULT_OFFSET_4_MAX_GAP_NUM), "UINT");
  opts.optopt("", "min_pow_of_2", &format!("A minimum power of 2 to calculate a gamma parameter (Uses {} by default)", DEFAULT_MIN_POW_OF_2), "FLOAT");
  opts.optopt("", "max_pow_of_2", &format!("A maximum power of 2 to calculate a gamma parameter (Uses {} by default)", DEFAULT_MAX_POW_OF_2), "FLOAT");
  opts.optopt("", "mix_weight", &format!("A mixture weight (Uses {} by default)", DEFAULT_MIX_WEIGHT), "FLOAT");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of all the threads of this computer by default)", "UINT");
  opts.optflag("b", "takes_bench", &format!("Compute for only gamma = {} to measure running time", GAMMA_4_BENCH));
  opts.optflag("q", "produces_access_probs", &format!("Also compute accessible probabilities"));
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
  let min_bpp = if matches.opt_present("min_base_pair_prob") {
    matches.opt_str("min_base_pair_prob").unwrap().parse().unwrap()
  } else {
    DEFAULT_MIN_BPP
  };
  let offset_4_max_gap_num = if matches.opt_present("offset_4_max_gap_num") {
    matches.opt_str("offset_4_max_gap_num").unwrap().parse().unwrap()
  } else {
    DEFAULT_OFFSET_4_MAX_GAP_NUM
  };
  let min_pow_of_2 = if matches.opt_present("min_pow_of_2") {
    matches.opt_str("min_pow_of_2").unwrap().parse().unwrap()
  } else {
    DEFAULT_MIN_POW_OF_2
  };
  let max_pow_of_2 = if matches.opt_present("max_pow_of_2") {
    matches.opt_str("max_pow_of_2").unwrap().parse().unwrap()
  } else {
    DEFAULT_MAX_POW_OF_2
  };
  let mix_weight = if matches.opt_present("mix_weight") {
    matches.opt_str("mix_weight").unwrap().parse().unwrap()
  } else {
    DEFAULT_MIX_WEIGHT
  };
  let takes_bench = matches.opt_present("b");
  let produces_access_probs = matches.opt_present("q");
  let outputs_probs = matches.opt_present("p");
  let num_of_threads = if matches.opt_present("t") {
    matches.opt_str("t").unwrap().parse().unwrap()
  } else {
    num_cpus::get() as NumOfThreads
  };
  let (cols, seq_ids) = read_sa_from_clustal_file(input_file_path);
  let sa_len = cols.len();
  let mut thread_pool = Pool::new(num_of_threads);
  if sa_len <= u8::MAX as usize {
    multi_threaded_consalifold::<u8>(&mut thread_pool, &cols, &seq_ids, offset_4_max_gap_num, min_bpp, produces_access_probs, output_dir_path, min_pow_of_2, max_pow_of_2, takes_bench, outputs_probs, mix_weight);
  } else {
    multi_threaded_consalifold::<u16>(&mut thread_pool, &cols, &seq_ids, offset_4_max_gap_num, min_bpp, produces_access_probs, output_dir_path, min_pow_of_2, max_pow_of_2, takes_bench, outputs_probs, mix_weight);
  }
}

fn multi_threaded_consalifold<T>(thread_pool: &mut Pool, cols: &Cols, seq_ids: &SeqIds, offset_4_max_gap_num: usize, min_bpp: Prob, produces_access_probs: bool, output_dir_path: &Path, min_pow_of_2: i32, max_pow_of_2: i32, takes_bench: bool, outputs_probs: bool, mix_weight: Prob)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Display + Sync + Send,
{
  let feature_score_sets = FeatureCountSets::load_trained_score_params();
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
      }
      if seq_lens[j] > 0 {
        sa.pos_map_sets[i][j] = T::from_usize(seq_lens[j]).unwrap();
      }
    }
  }
  for i in 0 .. num_of_rnas {
    fasta_records[i].seq.insert(0, PSEUDO_BASE);
    fasta_records[i].seq.push(PSEUDO_BASE);
    fasta_records[i].fasta_id = seq_ids[i].clone();
  }
  let mut rnaalifold_bpp_mat = SparseProbMat::<T>::default();
  let ref mut ref_2_rnaalifold_bpp_mat = rnaalifold_bpp_mat;
  let ref ref_2_sa = sa;
  let ref ref_2_fasta_records = fasta_records;
  let ref ref_2_feature_score_sets = feature_score_sets;
  thread_pool.scoped(|scope| {
    scope.execute(move || {
      *ref_2_rnaalifold_bpp_mat = rnaalifold_trained(ref_2_sa, ref_2_fasta_records, ref_2_feature_score_sets);
    });
  });
  let prob_mat_sets = consprob::<T>(thread_pool, &fasta_records, min_bpp, T::from_usize(offset_4_max_gap_num).unwrap(), produces_access_probs);
  let mix_bpp_mat = get_mix_bpp_mat(&prob_mat_sets, &rnaalifold_bpp_mat, &sa, mix_weight);
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  if takes_bench {
    let output_file_path = output_dir_path.join(&format!("gamma={}.sth", GAMMA_4_BENCH));
    compute_and_write_mea_css(&mix_bpp_mat, &sa, GAMMA_4_BENCH, &output_file_path, &fasta_records);
  } else {
    thread_pool.scoped(|scope| {
      for pow_of_2 in min_pow_of_2 .. max_pow_of_2 + 1 {
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
  if outputs_probs {
    write_prob_mat_sets::<T>(&output_dir_path, &prob_mat_sets, produces_access_probs);
  }
}

fn read_sa_from_clustal_file(clustal_file_path: &Path) -> (Cols, SeqIds) {
  let mut cols = Cols::new();
  let mut seq_ids = SeqIds::new();
  let reader_2_clustal_file = BufReader::new(File::open(clustal_file_path).unwrap());
  let mut seq_pointer = 0;
  let mut pos_pointer = 0;
  let mut are_seq_ids_read = false;
  for (i, string) in reader_2_clustal_file.lines().enumerate() {
    let string = string.unwrap();
    if i == 0 || string.len() == 0 || string.starts_with(" ") {
      if cols.len() > 0 {
        seq_pointer = 0;
        pos_pointer = cols.len();
        are_seq_ids_read = true;
      }
      continue;
    }
    let mut substrings = string.split_whitespace();
    let substring = substrings.next().unwrap();
    if !are_seq_ids_read {
      seq_ids.push(String::from(substring));
    }
    let substring = substrings.next().unwrap();
    if seq_pointer == 0 {
      for sa_char in substring.chars() {
        cols.push(vec![convert_char(sa_char as u8)]);
      }
      seq_pointer += 1;
    } else {
      for (j, sa_char) in substring.chars().enumerate() {
        cols[pos_pointer + j].push(convert_char(sa_char as u8));
      }
    }
  }
  (cols, seq_ids)
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
        // if sa.cols[i][k] == GAP || sa.cols[j][k] == GAP {continue;}
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
      let pos_pair = (T::from_usize(i).unwrap(), T::from_usize(j).unwrap());
      match rnaalifold_bpp_mat.get(&pos_pair) {
        Some(&rnaalifold_bpp) => {
          mix_bpp_mat[i][j] = if effective_num_of_rnas > 0 {
            mix_weight * mean_bpp / effective_num_of_rnas as Prob + (1. - mix_weight) * rnaalifold_bpp
          } else {
            rnaalifold_bpp
          };
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
