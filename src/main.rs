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
const DEFAULT_MIN_POW_OF_2: i32 = -4;
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
  opts.optopt("", "input_locarnap_pair_prob_matrix_file_path", "The path to an input file containing pairing probability matrices computed by the LocARNA-P algorithm", "STR");
  opts.optflag("u", "is_posterior_model", "Uses the posterior model to score secondary structures (Must specify the flag \"input_locarnap_pair_prob_matrix_file_path\")");
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
  let is_posterior_model = matches.opt_present("u");
  if is_posterior_model && !matches.opt_present("input_locarnap_pair_prob_matrix_file_path") {
    print_program_usage(&program_name, &opts);
    return;
  }
  let input_locarnap_bpp_mat_file_path = if matches.opt_present("input_locarnap_pair_prob_matrix_file_path") {
    matches.opt_str("input_locarnap_pair_prob_matrix_file_path").unwrap()
  } else {
    String::new()
  };
  let input_locarnap_bpp_mat_file_path = Path::new(&input_locarnap_bpp_mat_file_path);
  let takes_bench = matches.opt_present("b");
  let produces_access_probs = matches.opt_present("q") && !is_posterior_model;
  let outputs_probs = matches.opt_present("p");
  let num_of_threads = if matches.opt_present("t") {
    matches.opt_str("t").unwrap().parse().unwrap()
  } else {
    num_cpus::get() as NumOfThreads
  };
  let fasta_file_reader = Reader::from_file(Path::new(&input_file_path)).unwrap();
  let mut fasta_records = FastaRecords::new();
  let mut max_seq_len = 0;
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.unwrap();
    let mut seq = convert(fasta_record.seq());
    seq.insert(0, PSEUDO_BASE);
    seq.push(PSEUDO_BASE);
    let seq_len = seq.len();
    if seq_len > max_seq_len {
      max_seq_len = seq_len;
    }
    fasta_records.push(FastaRecord::new(String::from(fasta_record.id()), seq));
  }
  let mut thread_pool = Pool::new(num_of_threads);
  if max_seq_len <= u8::MAX as usize {
    multi_threaded_consalifold::<u8>(&mut thread_pool, &fasta_records, offset_4_max_gap_num, min_bpp, produces_access_probs, output_dir_path, min_pow_of_2, max_pow_of_2, takes_bench, outputs_probs, mix_weight, input_bpp_mat_file_path, input_locarnap_bpp_mat_file_path, is_posterior_model, input_sa_file_path);
  } else {
    multi_threaded_consalifold::<u16>(&mut thread_pool, &fasta_records, offset_4_max_gap_num, min_bpp, produces_access_probs, output_dir_path, min_pow_of_2, max_pow_of_2, takes_bench, outputs_probs, mix_weight, input_bpp_mat_file_path, input_locarnap_bpp_mat_file_path, is_posterior_model, input_sa_file_path);
  }
}

fn multi_threaded_consalifold<T>(thread_pool: &mut Pool, fasta_records: &FastaRecords, offset_4_max_gap_num: usize, min_bpp: Prob, produces_access_probs: bool, output_dir_path: &Path, min_pow_of_2: i32, max_pow_of_2: i32, takes_bench: bool, outputs_probs: bool, mix_weight: Prob, input_bpp_mat_file_path: &Path, input_locarnap_bpp_mat_file_path: &Path, is_posterior_model: bool, input_sa_file_path: &Path)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Display + Sync + Send,
{
  let (mut sa, seq_ids) = read_sa_from_clustal_file::<u16>(input_sa_file_path);
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
        sa.pos_map_sets[i][j] = seq_lens[j] as u16 - 1;
      }
    }
  }
  let sa_len = sa.cols.len();
  let prob_mat_sets = if is_posterior_model {
    read_locarnap_bpp_mats::<T>(&input_locarnap_bpp_mat_file_path, &seq_lens)
  } else {
    consprob::<T>(thread_pool, &fasta_records, min_bpp, T::from_usize(offset_4_max_gap_num).unwrap(), produces_access_probs)
  };
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
  let mix_bpp_mat = get_mix_bpp_mat::<T, u16>(&prob_mat_sets, &rnaalipfold_bpp_mat, &sa, mix_weight);
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  if takes_bench {
    let output_file_path = output_dir_path.join(&format!("gamma={}.sth", GAMMA_4_BENCH));
    compute_and_write_mea_css::<u16>(&mix_bpp_mat, &sa, GAMMA_4_BENCH, &output_file_path, &seq_ids);
  } else {
    thread_pool.scoped(|scope| {
      for pow_of_2 in min_pow_of_2 .. max_pow_of_2 + 1 {
        let gamma = (2. as Prob).powi(pow_of_2);
        let ref ref_2_mix_bpp_mat = mix_bpp_mat;
        let ref ref_2_sa = sa;
        let ref ref_2_seq_ids = seq_ids;
        let output_file_path = output_dir_path.join(&format!("gamma={}.sth", gamma));
        scope.execute(move || {
          compute_and_write_mea_css::<u16>(ref_2_mix_bpp_mat, ref_2_sa, gamma, &output_file_path, ref_2_seq_ids);
        });
      }
    });
  }
  if outputs_probs {
    write_prob_mat_sets::<T>(&output_dir_path, &prob_mat_sets, produces_access_probs);
  }
}

fn read_locarnap_bpp_mats<T>(input_locarnap_bpp_mat_file_path: &Path, seq_lens: &Vec<usize>) -> ProbMatSets<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Display + Sync + Send,
{
  let mut prob_mat_sets = Vec::new();
  let reader_2_locarnap_bpp_mat_file = BufReader::new(File::open(input_locarnap_bpp_mat_file_path).unwrap());
  for (i, string) in reader_2_locarnap_bpp_mat_file.split(b'>').enumerate() {
    let string = String::from_utf8(string.unwrap()).unwrap();
    if string.len() == 0 {continue;}
    let seq_len = seq_lens[i - 1] + 2;
    let mut prob_mats = PctStaProbMats::new(seq_len);
    match string.split('\n').skip(1).next() {
      Some(string) => {
        let string = string.trim();
        if string.len() > 0 {
          for substring in string.split(' ') {
            let subsubstrings = substring.split(',').map(|subsubstring| {String::from(subsubstring)}).collect::<Vec<String>>();
            let (m, n, bpp) = (subsubstrings[0].parse::<usize>().unwrap(), subsubstrings[1].parse::<usize>().unwrap(), subsubstrings[2].parse().unwrap());
            let pos_pair = (T::from_usize(m).unwrap(), T::from_usize(n).unwrap());
            prob_mats.bpp_mat.insert(pos_pair, bpp);
          }
        }
      }, None => {},
    }
    prob_mat_sets.push(prob_mats);
  }
  prob_mat_sets
}

fn read_sa_from_clustal_file<T>(clustal_file_path: &Path) -> (SeqAlign<T>, SeqIds)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Display + Sync + Send,
{
  let mut sa = SeqAlign::<T>::new();
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

fn get_mix_bpp_mat<T, U>(prob_mat_sets: &ProbMatSets<T>, rnaalipfold_bpp_mat: &ProbMat, sa: &SeqAlign<U>, mix_weight: Prob) -> ProbMat
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Display + Sync + Send,
  U: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Display + Sync + Send,
{
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
        let pos_pair = (T::from(sa.pos_map_sets[i][k]).unwrap() + T::one(), T::from(sa.pos_map_sets[j][k]).unwrap() + T::one());
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

fn compute_and_write_mea_css<T>(mix_bpp_mat: &ProbMat, sa: &SeqAlign<T>, gamma: Prob, output_file_path: &Path, seq_ids: &SeqIds)
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Display + Sync + Send,
{
  let mea_css = consalifold::<T>(mix_bpp_mat, gamma, sa);
  let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).unwrap());
  let mut buf_4_writer_2_output_file = format!("# STOCKHOLM 1.0\n");
  let sa_len = sa.cols.len();
  let max_seq_id_len = seq_ids.iter().map(|seq_id| {seq_id.len()}).max().unwrap();
  let num_of_rnas = sa.cols[0].len();
  for rna_id in 0 .. num_of_rnas {
    let ref seq_id = seq_ids[rna_id];
    buf_4_writer_2_output_file.push_str(seq_id);
    let mut stockholm_row = vec![' ' as Char; max_seq_id_len - seq_id.len() + 2];
    let mut sa_row = (0 .. sa_len).map(|x| {sa.cols[x][rna_id]}).collect::<Vec<Char>>();
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
