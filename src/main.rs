extern crate consalifold;
extern crate bio;
extern crate num_cpus;

use consalifold::*;
use std::env;
use std::path::Path;
use std::io::{BufReader, BufWriter};
use std::fs::File;
use bio::io::fasta::Reader;
use std::fs::create_dir;

type MeaCssStr = MeaSsStr;
type Strings = Vec<String>;

const DEFAULT_MIX_WEIGHT: Prob = 0.5;
const DEFAULT_MIN_POW_OF_2: i32 = -7;
const DEFAULT_MAX_POW_OF_2: i32 = 10;
const GAMMA_4_BENCH: Prob = 1.;

fn main() {
  let args = env::args().collect::<Vec<Arg>>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("i", "input_file_path", "The path to an input FASTA file that contains RNA sequences", "STR");
  opts.reqopt("a", "input_seq_align_file_path", "The path to an input CLUSTAL file containing a sequence alignment of RNA sequences", "STR");
  opts.reqopt("c", "input_pair_prob_matrix_file_path", "The path to an input file containing pairing probability matrices computed by the RNAalipfold algorithm", "STR");
  opts.reqopt("o", "output_dir_path", "The path to an output directory", "STR");
  opts.optopt("", "min_base_pair_prob", &format!("A minimum base-pairing-probability (Uses {} by default)", DEFAULT_MIN_BPP), "FLOAT");
  opts.optopt("", "offset_4_max_gap_num", &format!("An offset for maximum numbers of gaps (Uses {} by default)", DEFAULT_OFFSET_4_MAX_GAP_NUM), "UINT");
  opts.optopt("", "min_pow_of_2", &format!("A minimum power of 2 to calculate a gamma parameter (Uses {} by default)", DEFAULT_MIN_POW_OF_2), "FLOAT");
  opts.optopt("", "max_pow_of_2", &format!("A maximum power of 2 to calculate a gamma parameter (Uses {} by default)", DEFAULT_MAX_POW_OF_2), "FLOAT");
  opts.optopt("", "mix_weight", &format!("A mixture weight (Uses {} by default)", DEFAULT_MIX_WEIGHT), "FLOAT");
  opts.optflag("u", "uses_bpp_score", "Uses base-pairing probabilities as scores of secondary structures (Not recommended due to poor accuracy)");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of all the threads of this computer by default)", "UINT");
  opts.optflag("b", "takes_bench", &format!("Compute for only gamma = {} to measure running time", GAMMA_4_BENCH));
  opts.optflag("a", "produces_access_probs", &format!("Also compute accessible probabilities"));
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
  let input_sa_file_path = matches.opt_str("a").unwrap();
  let input_sa_file_path = Path::new(&input_sa_file_path);
  let input_bpp_mat_file_path = matches.opt_str("c").unwrap();
  let input_bpp_mat_file_path = Path::new(&input_bpp_mat_file_path);
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
  let uses_bpps = matches.opt_present("u");
  let takes_bench = matches.opt_present("b");
  let produces_access_probs = matches.opt_present("a");
  let outputs_probs = matches.opt_present("p");
  let num_of_threads = if matches.opt_present("t") {
    matches.opt_str("t").unwrap().parse().unwrap()
  } else {
    num_cpus::get() as NumOfThreads
  };
  let fasta_file_reader = Reader::from_file(Path::new(&input_file_path)).unwrap();
  let mut fasta_records = FastaRecords::new();
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.unwrap();
    let mut seq = convert(fasta_record.seq());
    seq.insert(0, PSEUDO_BASE);
    seq.push(PSEUDO_BASE);
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let mut thread_pool = Pool::new(num_of_threads);
  let prob_mat_sets = consprob(&mut thread_pool, &fasta_records, min_bpp, offset_4_max_gap_num, uses_bpps, produces_access_probs);
  let (mut sa, seq_ids) = read_sa_from_clustal_file(input_sa_file_path);
  let num_of_rnas = sa.cols[0].len();
  let mut seq_lens = vec![0 as usize; num_of_rnas];
  let num_of_cols = sa.cols.len();
  sa.pos_map_sets = vec![vec![0; num_of_rnas]; num_of_cols];
  for i in 0 .. num_of_cols {
    for j in 0 .. num_of_rnas {
      let base = sa.cols[i][j];
      if base != GAP {
        seq_lens[j] += 1;
      }
      if seq_lens[j] > 0 {
        sa.pos_map_sets[i][j] = seq_lens[j] as Pos - 1;
      }
    }
  }
  let sa_len = sa.cols.len();
  let mut rnaalipfold_bpp_mat = vec![vec![0.; sa_len]; sa_len];
  let reader_2_input_bpp_mat_file = BufReader::new(File::open(input_bpp_mat_file_path).unwrap());
  for (i, vec) in reader_2_input_bpp_mat_file.split(b'\n').enumerate() {
    let vec = vec.unwrap();
    let substrings = unsafe {String::from_utf8_unchecked(vec).split_whitespace().map(|string| {String::from(string)}).collect::<Strings>()};
    for subsubstring in &substrings[2 ..] {
      let subsubsubstrings = subsubstring.split(":").collect::<Vec<&str>>();
      let j = subsubsubstrings[0].parse::<usize>().unwrap() - 1;
      rnaalipfold_bpp_mat[i][j] = subsubsubstrings[1].parse().unwrap();
    }
  }
  let mix_bpp_mat = get_mix_bpp_mat(&prob_mat_sets, &rnaalipfold_bpp_mat, &sa, mix_weight);
  let mix_upp_mat = get_mix_upp_mat(&prob_mat_sets, &rnaalipfold_bpp_mat, &sa, mix_weight);
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  if takes_bench {
    let output_file_path = output_dir_path.join(&format!("gamma={}.sth", GAMMA_4_BENCH));
    compute_and_write_mea_css(&mix_bpp_mat, &mix_upp_mat, &sa, GAMMA_4_BENCH, &output_file_path, &seq_ids);
  } else {
    thread_pool.scoped(|scope| {
      for pow_of_2 in min_pow_of_2 .. max_pow_of_2 + 1 {
        let gamma = (2. as Prob).powi(pow_of_2);
        let ref ref_2_mix_bpp_mat = mix_bpp_mat;
        let ref ref_2_mix_upp_mat = mix_upp_mat;
        let ref ref_2_sa = sa;
        let ref ref_2_seq_ids = seq_ids;
        let output_file_path = output_dir_path.join(&format!("gamma={}.sth", gamma));
        scope.execute(move || {
          compute_and_write_mea_css(ref_2_mix_bpp_mat, ref_2_mix_upp_mat, ref_2_sa, gamma, &output_file_path, ref_2_seq_ids);
        });
      }
    });
  }
  if outputs_probs {
    write_prob_mat_sets(&output_dir_path, &prob_mat_sets, produces_access_probs);
  }
}

#[inline]
fn read_sa_from_clustal_file(clustal_file_path: &Path) -> (SeqAlign, SeqIds) {
  let mut sa = SeqAlign::new();
  let mut seq_ids = SeqIds::new();
  let reader_2_clustal_file = BufReader::new(File::open(clustal_file_path).unwrap());
  let mut seq_pointer = 0;
  let mut pos_pointer = 0;
  let mut are_seq_ids_read = false;
  for (i, string) in reader_2_clustal_file.lines().enumerate() {
    let string = string.unwrap();
    if i == 0 || string.len() == 0 || string.starts_with(" ") {
      if sa.cols.len() > 0 {
        seq_pointer = 0;
        pos_pointer = sa.cols.len();
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
        sa.cols.push(vec![sa_char as Char]);
      }
      seq_pointer += 1;
    } else {
      for (j, sa_char) in substring.chars().enumerate() {
        sa.cols[pos_pointer + j].push(sa_char as Char);
      }
    }
  }
  (sa, seq_ids)
}

#[inline]
fn get_mix_bpp_mat(prob_mat_sets: &ProbMatSets, rnaalipfold_bpp_mat: &ProbMat, sa: &SeqAlign, mix_weight: Prob) -> ProbMat {
  let sa_len = sa.cols.len();
  let num_of_rnas = sa.cols[0].len();
  let mut mix_bpp_mat = vec![vec![0.; sa_len]; sa_len];
  for i in 0 .. sa_len {
    for j in i + 1 .. sa_len {
      let mut mean_bpp = 0.;
      let mut effective_num_of_rnas = 0;
      for k in 0 .. num_of_rnas {
        if sa.cols[i][k] == GAP || sa.cols[j][k] == GAP {continue;}
        let ref bpp_mat = prob_mat_sets[k].bpp_mat;
        let pos_pair = (sa.pos_map_sets[i][k] + 1, sa.pos_map_sets[j][k] + 1);
        match bpp_mat.get(&pos_pair) {
          Some(&bpp) => {
            mean_bpp += bpp;
            effective_num_of_rnas += 1;
          }, None => {},
        }
      }
      mix_bpp_mat[i][j] = if effective_num_of_rnas > 0 {
        mix_weight * mean_bpp / effective_num_of_rnas as Prob + (1. - mix_weight) * rnaalipfold_bpp_mat[i][j]
      } else {
        rnaalipfold_bpp_mat[i][j]
      };
    }
  }
  mix_bpp_mat
}

#[inline]
fn get_mix_upp_mat(prob_mat_sets: &ProbMatSets, rnaalipfold_bpp_mat: &ProbMat, sa: &SeqAlign, mix_weight: Prob) -> Probs {
  let sa_len = sa.cols.len();
  let num_of_rnas = sa.cols[0].len();
  let mut mix_upp_mat = vec![0.; sa_len];
  for i in 0 .. sa_len {
    let mut rnaalipfold_upp = 1.;
    for j in 0 .. sa_len {
      if i == j {continue;}
      let pos_pair = if i < j {(i, j)} else {(j, i)};
      rnaalipfold_upp -= rnaalipfold_bpp_mat[pos_pair.0][pos_pair.1];
    }
    let mut mean_upp = 0.;
    let mut effective_num_of_rnas = 0;
    for j in 0 .. num_of_rnas {
      if sa.cols[i][j] == GAP {continue;}
      mean_upp += prob_mat_sets[j].upp_mat[sa.pos_map_sets[i][j] as usize + 1];
      effective_num_of_rnas += 1;
    }
    mix_upp_mat[i] = if effective_num_of_rnas > 0 {
      mix_weight * mean_upp / num_of_rnas as Prob + (1. - mix_weight) * rnaalipfold_upp
    } else {
      rnaalipfold_upp
    };
  }
  mix_upp_mat
}

fn compute_and_write_mea_css(mix_bpp_mat: &ProbMat, mix_upp_mat: &Probs, sa: &SeqAlign, gamma: Prob, output_file_path: &Path, seq_ids: &SeqIds) {
  let mea_css = consalifold(mix_bpp_mat, mix_upp_mat, gamma, sa);
  let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).unwrap());
  let mut buf_4_writer_2_output_file = format!("# STOCKHOLM 1.0\n");
  let sa_len = sa.cols.len();
  let max_seq_id_len = seq_ids.iter().map(|seq_id| {seq_id.len()}).max().unwrap();
  let num_of_rnas = sa.cols[0].len();
  for rna_id in 0 .. num_of_rnas {
    let ref seq_id = seq_ids[rna_id];
    buf_4_writer_2_output_file.push_str(seq_id);
    let mut stockholm_row = vec![' ' as Char; max_seq_id_len - seq_id.len() + 2];
    let mut sa_row = (0 .. sa_len).map(|x| {sa.cols[x][rna_id]}).collect::<Vec<u8>>();
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

#[inline]
fn get_mea_css_str(mea_css: &MeaCss, sa_len: usize) -> MeaCssStr {
  let mut mea_css_str = vec![UNPAIRING_BASE; sa_len];
  for &(i, j) in &mea_css.bpa_pos_pairs {
    mea_css_str[i as usize] = BASE_PAIRING_LEFT_BASE;
    mea_css_str[j as usize] = BASE_PAIRING_RIGHT_BASE;
  }
  mea_css_str
}
