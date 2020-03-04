extern crate phyloalifold;
extern crate getopts;

use phyloalifold::*;
use getopts::Options;
use std::env;
use std::path::Path;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::fs::File;

type Arg = String;
type Strings = Vec<String>;
type MeaCssStr = MeaSsStr;

const DEFAULT_GAMMA: Prob = 1.;
const DEFAULT_PROB_WEIGHT: Prob = 0.5;
const VERSION: &'static str = "0.1.0";

fn main() {
  let args = env::args().collect::<Vec<Arg>>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("a", "input_seq_align_file_path", "The path to an input CLUSTAL file containing a sequence alignment of RNA sequences", "STR");
  opts.reqopt("p", "input_base_pair_prob_matrix_file_path", "The path to an input file containing base-pairing probability matrices", "STR");
  opts.reqopt("q", "input_unpair_prob_matrix_file_path", "The path to an input file containing unpairing probability matrices", "STR");
  opts.reqopt("r", "input_base_pair_prob_matrix_file_path_2", "The path to an input file containing base-pairing probability matrices computed by the CentroidAlifold algorithm", "STR");
  opts.reqopt("o", "output_file_path", "The path to an output STOCKHOLM file which will contain estimated consensus secondary structures", "STR");
  opts.optopt("", "gamma", &format!("An MEA gamma (Uses {} by default)", DEFAULT_GAMMA), "FLOAT");
  opts.optopt("", "prob_weight", &format!("A probability weight (Uses {} by default)", DEFAULT_PROB_WEIGHT), "FLOAT");
  opts.optflag("h", "help", "Print a help menu");
  let opts = match opts.parse(&args[1 ..]) {
    Ok(opt) => {opt}
    Err(failure) => {print_program_usage(&program_name, &opts); panic!(failure.to_string())}
  };
  let input_sa_file_path = opts.opt_str("a").unwrap();
  let input_sa_file_path = Path::new(&input_sa_file_path);
  let input_bpp_mat_file_path = opts.opt_str("p").unwrap();
  let input_bpp_mat_file_path = Path::new(&input_bpp_mat_file_path);
  let input_upp_mat_file_path = opts.opt_str("q").unwrap();
  let input_upp_mat_file_path = Path::new(&input_upp_mat_file_path);
  let input_bpp_mat_file_path_2 = opts.opt_str("r").unwrap();
  let input_bpp_mat_file_path_2 = Path::new(&input_bpp_mat_file_path_2);
  let output_file_path = opts.opt_str("o").unwrap();
  let output_file_path = Path::new(&output_file_path);
  let gamma = if opts.opt_present("gamma") {
    opts.opt_str("gamma").unwrap().parse().unwrap()
  } else {
    DEFAULT_GAMMA
  };
  let prob_weight = if opts.opt_present("prob_weight") {
    opts.opt_str("prob_weight").unwrap().parse().unwrap()
  } else {
    DEFAULT_PROB_WEIGHT
  };
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
  let mut bpp_mats_with_rna_ids = vec![SparseProbMat::default(); num_of_rnas];
  let mut reader_2_input_bpp_mat_file = BufReader::new(File::open(input_bpp_mat_file_path).unwrap());
  let mut buf_4_reader_2_input_bpp_mat_file = Vec::new();
  for _ in 0 .. 2 {
    let _ = reader_2_input_bpp_mat_file.read_until(b'>', &mut buf_4_reader_2_input_bpp_mat_file);
  }
  for (i, vec) in reader_2_input_bpp_mat_file.split(b'>').enumerate() {
    if i == num_of_rnas {break;}
    let vec = vec.unwrap();
    let substrings = unsafe {String::from_utf8_unchecked(vec).split_whitespace().map(|string| {String::from(string)}).collect::<Strings>()};
    let rna_id = substrings[0].parse::<RnaId>().unwrap();
    let ref mut bpp_mat = bpp_mats_with_rna_ids[rna_id];
    *bpp_mat = SparseProbMat::default();
    for subsubstring in &substrings[1 ..] {
      let subsubsubstrings = subsubstring.split(",").collect::<Vec<&str>>();
      let pos_pair = (
        subsubsubstrings[0].parse::<Pos>().unwrap(),
        subsubsubstrings[1].parse::<Pos>().unwrap(),
        );
      bpp_mat.insert(pos_pair, subsubsubstrings[2].parse().unwrap());
    }
  }
  let mean_bpp_mat = get_mean_bpp_mat(&bpp_mats_with_rna_ids, &sa);
  let mut upp_mats_with_rna_ids = vec![Probs::new(); num_of_rnas];
  let mut reader_2_input_upp_mat_file = BufReader::new(File::open(input_upp_mat_file_path).unwrap());
  let mut buf_4_reader_2_input_upp_mat_file = Vec::new();
  for _ in 0 .. 2 {
    let _ = reader_2_input_upp_mat_file.read_until(b'>', &mut buf_4_reader_2_input_upp_mat_file);
  }
  for (i, vec) in reader_2_input_upp_mat_file.split(b'>').enumerate() {
    if i == num_of_rnas {break;}
    let vec = vec.unwrap();
    let substrings = unsafe {String::from_utf8_unchecked(vec).split_whitespace().map(|string| {String::from(string)}).collect::<Strings>()};
    let rna_id = substrings[0].parse::<RnaId>().unwrap();
    let seq_len = seq_lens[rna_id];
    let ref mut upp_mat = upp_mats_with_rna_ids[rna_id];
    *upp_mat = vec![0.; seq_len];
    for subsubstring in &substrings[1 ..] {
      let subsubsubstrings = subsubstring.split(",").collect::<Vec<&str>>();
      upp_mat[subsubsubstrings[0].parse::<usize>().unwrap()] = subsubsubstrings[1].parse().unwrap();
    }
  }
  let mean_upp_mat = get_mean_upp_mat(&upp_mats_with_rna_ids, &sa);
  let sa_len = sa.cols.len();
  let num_of_rnas = sa.cols[0].len();
  let mut centroidalifold_bpp_mat = vec![vec![0.; sa_len]; sa_len];
  let reader_2_input_bpp_mat_file_2 = BufReader::new(File::open(input_bpp_mat_file_path_2).unwrap());
  for (i, vec) in reader_2_input_bpp_mat_file_2.split(b'\n').enumerate() {
    let vec = vec.unwrap();
    let substrings = unsafe {String::from_utf8_unchecked(vec).split_whitespace().map(|string| {String::from(string)}).collect::<Strings>()};
    for subsubstring in &substrings[2 ..] {
      let subsubsubstrings = subsubstring.split(":").collect::<Vec<&str>>();
      let j = subsubsubstrings[0].parse::<usize>().unwrap() - 1;
      centroidalifold_bpp_mat[i][j] = subsubsubstrings[1].parse().unwrap();
    }
  }
  let mea_css = neoalifold(&mean_bpp_mat, &mean_upp_mat, &centroidalifold_bpp_mat, gamma, prob_weight, &sa);
  let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).unwrap());
  let mut buf_4_writer_2_output_file = format!("# STOCKHOLM 1.0\n#=GF CC The version {} of the NeoAliFold program\n#=GF CC The path to the input CLUSTAL file for computing the consensus secondary structure (= CSS) in this file = \"{}\"\n#=GF CC The path to the input base-pairing probability matrix file for computing this structure = \"{}\"\n#=GF CC The path to the input unpairing probability matrix file for computing this structure = \"{}\"\n#=GF CC The values of the parameters used for computing this structure are as follows\n#=GF CC \"gamma\" = {}\n", VERSION, input_sa_file_path.display(), input_bpp_mat_file_path.display(), input_upp_mat_file_path.display(), gamma);
  let sa_len = sa.cols.len();
  let max_seq_id_len = seq_ids.iter().map(|seq_id| {seq_id.len()}).max().unwrap();
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
fn print_program_usage(program_name: &str, opts: &Options) {
  let program_usage = format!("The usage of this program: {} [options]", program_name);
  print!("{}", opts.usage(&program_usage));
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
fn get_mean_bpp_mat(bpp_mats_with_rna_ids: &ProbMatsWithRnaIds, sa: &SeqAlign) -> ProbMat {
  let sa_len = sa.cols.len();
  let num_of_rnas = sa.cols[0].len();
  let mut mean_bpp_mat = vec![vec![0.; sa_len]; sa_len];
  for i in 0 .. sa_len {
    for j in i + 1 .. sa_len {
      let mut mean_bpp = 0.;
      let mut effective_num_of_rnas = 0;
      for k in 0 .. num_of_rnas {
        if sa.cols[i][k] == GAP || sa.cols[j][k] == GAP {continue;}
        let ref bpp_mat = bpp_mats_with_rna_ids[k];
        let pos_pair = (sa.pos_map_sets[i][k], sa.pos_map_sets[j][k]);
        if bpp_mat.contains_key(&pos_pair) {
          mean_bpp += bpp_mat[&pos_pair];
          effective_num_of_rnas += 1;
        }
      }
      if effective_num_of_rnas > 0 {
        mean_bpp_mat[i][j] = mean_bpp / effective_num_of_rnas as Prob;
      }
    }
  }
  mean_bpp_mat
}

#[inline]
fn get_mean_upp_mat(upp_mats_with_rna_ids: &ProbsWithRnaIds, sa: &SeqAlign) -> Probs {
  let sa_len = sa.cols.len();
  let num_of_rnas = sa.cols[0].len();
  let mut mean_upp_mat = vec![0.; sa_len];
  for i in 0 .. sa_len {
    let mut mean_upp = 0.;
    let mut effective_num_of_rnas = 0;
    for j in 0 .. num_of_rnas {
      if sa.cols[i][j] == GAP {continue;}
      mean_upp += upp_mats_with_rna_ids[j][sa.pos_map_sets[i][j] as usize];
      effective_num_of_rnas += 1;
    }
    if effective_num_of_rnas > 0 {
      mean_upp_mat[i] = mean_upp / num_of_rnas as Prob;
    }
  }
  mean_upp_mat
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
