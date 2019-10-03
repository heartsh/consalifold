extern crate neoalifold;
extern crate getopts;
/* extern crate scoped_threadpool;
extern crate bio;
extern crate num_cpus;
extern crate petgraph;
extern crate rand; */

use neoalifold::*;
use getopts::Options;
// use self::scoped_threadpool::Pool;
use std::env;
use std::path::Path;
// use bio::io::fasta::Reader;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::fs::File;
/* use petgraph::{Outgoing, Incoming};
use rand::Rng; */

type Arg = String;
// type NumOfThreads = u32;
type Strings = Vec<String>;
// type MeaCsss = Vec<MeaCss>;
type MeaCssStr = MeaSsStr;

const DEFAULT_GAMMA: Prob = 1.;
// const DEFAULT_NUM_OF_ITERATIVE_REFINEMENTS: usize = 5;
const VERSION: &'static str = "0.1.0";
const PSEUDO_BASE: Char = '$' as Char;

fn main() {
  let args = env::args().collect::<Vec<Arg>>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("a", "input_seq_align_file_path", "The path to an input CLUSTAL file containing a sequence alignment of RNA sequences", "STR");
  opts.reqopt("p", "input_base_pair_align_prob_matrix_file_path", "The path to an input file containing base pair alignment probability matrices", "STR");
  opts.reqopt("q", "input_unpair_prob_matrix_file_path", "The path to an input file containing unpairing probability matrices", "STR");
  opts.reqopt("o", "output_file_path", "The path to an output STOCKHOLM file which will contain estimated consensus secondary structures", "STR");
  opts.optopt("", "gamma", &format!("An MEA gamma (Uses {} by default)", DEFAULT_GAMMA), "FLOAT");
  // opts.optopt("", "num_of_iterative_refinements", &format!("The number of iterative refinements (Uses {} by default)", DEFAULT_NUM_OF_ITERATIVE_REFINEMENTS), "UINT");
  // opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of all the threads of this computer by default)", "UINT");
  opts.optflag("h", "help", "Print a help menu");
  let opts = match opts.parse(&args[1 ..]) {
    Ok(opt) => {opt}
    Err(failure) => {print_program_usage(&program_name, &opts); panic!(failure.to_string())}
  };
  let input_sa_file_path = opts.opt_str("a").expect("Failed to get the path to an input CLUSTAL file containing a sequence alignment of RNA sequences from command arguments.");
  let input_sa_file_path = Path::new(&input_sa_file_path);
  let input_bpap_mat_file_path = opts.opt_str("p").expect("Failed to get the path to an input file containing base pair alignment probability matrices from command arguments.");
  let input_bpap_mat_file_path = Path::new(&input_bpap_mat_file_path);
  let input_upp_mat_file_path = opts.opt_str("q").expect("Failed to get the path to an input file containing unpairing probability matrices from command arguments.");
  let input_upp_mat_file_path = Path::new(&input_upp_mat_file_path);
  let output_file_path = opts.opt_str("o").expect("Failed to get the path to an output file which will contain estimated consensus secondary structures from command arguments.");
  let output_file_path = Path::new(&output_file_path);
  let gamma = if opts.opt_present("gamma") {
    opts.opt_str("gamma").expect("Failed to get an MEA gamma from command arguments.").parse().expect("Failed to parse an MEA gamma.")
  } else {
    DEFAULT_GAMMA
  };
  /* let num_of_iterative_refinements = if opts.opt_present("num_of_iterative_refinements") {
    opts.opt_str("num_of_iterative_refinements").expect("Failed to get the number of iterative refinements from command arguments.").parse().expect("Failed to parse the number of iterative refinements.")
  } else {
    DEFAULT_NUM_OF_ITERATIVE_REFINEMENTS
  }; */
  /* let num_of_threads = if opts.opt_present("t") {
    opts.opt_str("t").expect("Failed to get the number of threads in multithreading from command arguments.").parse().expect("Failed to parse the number of threads in multithreading.")
  } else {
    num_cpus::get() as NumOfThreads
  }; */
  // let fasta_file_reader = Reader::from_file(Path::new(&input_fasta_file_path)).expect("Failed to set a FASTA file reader.");
  let (mut sa, seq_ids) = read_sa_from_clustal_file(input_sa_file_path);
  let num_of_rnas = sa.cols[0].len();
  let mut seq_lens = vec![0; num_of_rnas];
  let num_of_cols = sa.cols.len();
  sa.pos_map_sets = vec![vec![0; num_of_rnas]; num_of_cols];
  for i in 0 .. num_of_cols {
    for j in 0 .. num_of_rnas {
      let base = sa.cols[i][j];
      if base != GAP {
        seq_lens[j] += 1;
      }
      if seq_lens[j] > 0 {
        sa.pos_map_sets[i][j] = seq_lens[j] - 1;
      }
    }
  }
  /* for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.expect("Failed to read a FASTA record.");
    let seq = unsafe {from_utf8_unchecked(fasta_record.seq()).to_uppercase().as_bytes().iter().filter(|&&base| {is_rna_base(base)}).map(|&base| {base}).collect::<Seq>()};
    let seq_len = seq.len();
    fasta_records.push((String::from(fasta_record.id()), seq, seq_len));
  } */
  // let num_of_fasta_records = fasta_records.len();
  let mut bpap_mats_with_rna_id_pairs = Prob4dMatsWithRnaIdPairs::default();
  for rna_id_1 in 0 .. num_of_rnas {
    for rna_id_2 in rna_id_1 + 1 .. num_of_rnas {
      bpap_mats_with_rna_id_pairs.insert((rna_id_1, rna_id_2), Prob4dMat::default());
    }
  }
  let mut reader_2_input_bpap_mat_file = BufReader::new(File::open(input_bpap_mat_file_path).expect("Failed to read an input file."));
  let mut buf_4_reader_2_input_bpap_mat_file = Vec::new();
  for _ in 0 .. 2 {
    let _ = reader_2_input_bpap_mat_file.read_until(b'>', &mut buf_4_reader_2_input_bpap_mat_file);
  }
  for (i, vec) in reader_2_input_bpap_mat_file.split(b'>').enumerate() {
    if i == num_of_rnas * (num_of_rnas - 1) / 2 {break;}
    let vec = vec.expect("Failed to read an input file.");
    let substrings = unsafe {String::from_utf8_unchecked(vec).split_whitespace().map(|string| {String::from(string)}).collect::<Strings>()};
    let rna_id_pair = substrings[0].split(',').map(|string| {String::from(string)}).collect::<Strings>();
    let rna_id_pair = (
      rna_id_pair[0].parse::<RnaId>().expect("Failed to parse an RNA ID."),
      rna_id_pair[1].parse::<RnaId>().expect("Failed to parse an RNA ID."),
    );
    let seq_len_pair = (seq_lens[rna_id_pair.0], seq_lens[rna_id_pair.1]);
    let bpap_mat = bpap_mats_with_rna_id_pairs.get_mut(&rna_id_pair).expect("Failed to get an element of a hash map.");
    bpap_mat.insert((0, seq_len_pair.0 + 1, 0, seq_len_pair.1 + 1), 1.);
    for subsubstring in &substrings[1 ..] {
      let subsubsubstrings = subsubstring.split(",").collect::<Vec<&str>>();
      bpap_mat.insert((
        subsubsubstrings[0].parse::<Pos>().expect("Failed to parse an index.") + 1,
        subsubsubstrings[1].parse::<Pos>().expect("Failed to parse an index.") + 1,
        subsubsubstrings[2].parse::<Pos>().expect("Failed to parse an index.") + 1,
        subsubsubstrings[3].parse::<Pos>().expect("Failed to parse an index.") + 1),
        subsubsubstrings[4].parse().expect("Failed to parse a base pair alignment probability."),
      );
    }
  }
  let mut upp_mats_with_rna_ids = vec![Probs::new(); num_of_rnas];
  let mut reader_2_input_upp_mat_file = BufReader::new(File::open(input_upp_mat_file_path).expect("Failed to read an input file."));
  let mut buf_4_reader_2_input_upp_mat_file = Vec::new();
  for _ in 0 .. 2 {
    let _ = reader_2_input_upp_mat_file.read_until(b'>', &mut buf_4_reader_2_input_upp_mat_file);
  }
  for (i, vec) in reader_2_input_upp_mat_file.split(b'>').enumerate() {
    if i == num_of_rnas {break;}
    let vec = vec.expect("Failed to read an input file.");
    let substrings = unsafe {String::from_utf8_unchecked(vec).split_whitespace().map(|string| {String::from(string)}).collect::<Strings>()};
    let rna_id = substrings[0].parse::<RnaId>().expect("Failed to parse an RNA ID.");
    let seq_len = seq_lens[rna_id];
    let ref mut upp_mat = upp_mats_with_rna_ids[rna_id];
    *upp_mat = vec![0.; seq_len + 2];
    for subsubstring in &substrings[1 ..] {
      let subsubsubstrings = subsubstring.split(",").collect::<Vec<&str>>();
      upp_mat[subsubsubstrings[0].parse::<Pos>().expect("Failed to parse an index.") + 1] = subsubsubstrings[1].parse().expect("Failed to parse a base pair alignment probability.");
    }
  }
  let mean_upp_mat = get_mean_upp_mat(&upp_mats_with_rna_ids, &sa);
  let mea_css = neoalifold(&bpap_mats_with_rna_id_pairs, &mean_upp_mat, gamma, &sa);
  let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).expect("Failed to create an output file."));
  // let mut buf_4_writer_2_output_file = format!("; The version {} of the NeoAliFold program.\n; The path to the input CLUSTAL file for computing the consensus secondary structure (= CSS) in this file = \"{}\".\n; The path to the input base pair alignment probability matrix file for computing this structure = \"{}\".\n; The path to the input unpairing probability matrix file for computing this structure = \"{}\".\n; The values of the parameters used for computing this structure are as follows.\n; \"gamma\" = {}.\n\n", VERSION, input_sa_file_path.display(), input_bpap_mat_file_path.display(), input_upp_mat_file_path.display(), gamma);
  let mut buf_4_writer_2_output_file = format!("# STOCKHOLM 1.0\n#=GF CC The version {} of the NeoAliFold program\n#=GF CC The path to the input CLUSTAL file for computing the consensus secondary structure (= CSS) in this file = \"{}\"\n#=GF CC The path to the input base pair alignment probability matrix file for computing this structure = \"{}\"\n#=GF CC The path to the input unpairing probability matrix file for computing this structure = \"{}\"\n#=GF CC The values of the parameters used for computing this structure are as follows\n#=GF CC \"gamma\" = {}\n", VERSION, input_sa_file_path.display(), input_bpap_mat_file_path.display(), input_upp_mat_file_path.display(), gamma);
  let sa_len = sa.cols.len();
  let max_seq_id_len = seq_ids.iter().map(|seq_id| {seq_id.len()}).max().expect("Failed to find the maximum length of sequence IDs.");
  for rna_id in 0 .. num_of_rnas {
    let ref seq_id = seq_ids[rna_id];
    buf_4_writer_2_output_file.push_str(seq_id);
    let mut stockholm_row = vec![' ' as Char; max_seq_id_len - seq_id.len() + 2];
    let mut sa_row = (1 .. sa_len - 1).map(|x| {sa.cols[x][rna_id]}).collect::<Seq>();
    stockholm_row.append(&mut sa_row);
    let mut stockholm_row = unsafe {from_utf8_unchecked(&stockholm_row)};
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
  /* let mut mea_csss = vec![MeaCss::new(); num_of_fasta_records];
  for i in 0 .. num_of_fasta_records {
    let ref mut mea_css = mea_csss[i];
    mea_css.seq_num = 1;
    let ref seq = fasta_records[i].1;
    let seq_len = fasta_records[i].2;
    mea_css.col_num = seq_len + 2;
    for j in 0 .. seq_len + 2 {
      // mea_css.rna_id_pos_pair_seqs_with_cols.insert(j, vec![(i, j)]);
      mea_css.pos_seqs_with_cols.push(vec![j as FloatPos]);
      /* for k in j + MIN_HL_LEN + 1 .. seq_len + 2 {
        let pos_pair = (j, k);
        if !(j == 0 && k == seq_len + 1) && !(j > 0 && k < seq_len + 1 && is_canonical(&(seq[j - 1], seq[k - 1]))) {continue;}
        // mea_css.rna_id_pos_triple_seqs_with_col_pairs.insert(pos_pair, vec![(i, pos_pair.0, pos_pair.1)]);
        mea_css.pos_pair_seqs.insert(pos_pair, vec![(i, pos_pair.0, pos_pair.1)]);
      } */
    }
    mea_css.rna_ids = vec![i];
  }
  let mut mea_mat = SparseMeaMat::default();
  for rna_id_1 in 0 .. num_of_fasta_records {
    for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
      mea_mat.insert((rna_id_1, rna_id_2), NEG_INFINITY);
    }
  }
  let mut thread_pool = Pool::new(num_of_threads);
  thread_pool.scoped(|scope| {
    for (rna_id_pair, mea_mat) in mea_mat.iter_mut() {
      let mea_css_pair = (&mea_csss[rna_id_pair.0], &mea_csss[rna_id_pair.1]);
      let ref ref_2_bpap_mats_with_rna_id_pairs = bpap_mats_with_rna_id_pairs;
      let ref ref_2_upp_mats_with_rna_ids = upp_mats_with_rna_ids;
      scope.execute(move || {
        *mea_mat = get_mea_consensus_ss(&mea_css_pair, gamma_plus_1, ref_2_bpap_mats_with_rna_id_pairs, ref_2_upp_mats_with_rna_ids).mea;
      });
    }
  }); */
  /* println!("Computed MEA matrix.");
  let guide_tree = get_guide_tree(&mea_mat, num_of_fasta_records);
  let root_node = guide_tree.externals(Incoming).next().expect("Failed to get the root node of a guide tree.");
  let mea_css = get_mea_css_of_node(&guide_tree, &root_node, &mea_csss, gamma_plus_1, &bpap_mats_with_rna_id_pairs, &upp_mats_with_rna_ids, num_of_fasta_records, num_of_iterative_refinements);
  let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).expect("Failed to create an output file."));
  let mut buf_4_writer_2_output_file = format!("; The version {} of the MEARCOF program.\n; The path to the input FASTA file for computing the consensus secondary structures (= CSSs) in this file = \"{}\".\n; The path to the input base pair alignment matrix file for computing these structures = \"{}\".\n; The values of the parameters used for computing these structures are as follows.\n; \"gamma\" = {}, \"num_of_iterative_refinements\" = {}, \"num_of_threads\" = {}.\n; Each row is with pairs of the ID of each of RNA sequences and corresponding position pairs.\n\nmea = {}\n", VERSION, input_fasta_file_path.display(), input_bpap_mat_file_path.display(), gamma_plus_1 - 1., num_of_iterative_refinements, num_of_threads, mea_css.mea); */
  // let seq_num = mea_css.seq_num;
  /* for (col_pair, rna_id_pos_triples) in mea_css.rna_id_pos_triple_seqs_with_col_pairs.iter() {
    if col_pair.0 == 0 {continue;}
    let mut buf_4_col_pair = String::new();
    let num_of_rna_id_pos_triples = rna_id_pos_triples.len();
    for (i, &(rna_id, pos_1, pos_2)) in rna_id_pos_triples.iter().enumerate() {
      buf_4_col_pair.push_str(&if i < num_of_rna_id_pos_triples - 1 {
        format!("{}:{},{} ", rna_id, pos_1 - 1, pos_2 - 1)
      } else {
        format!("{}:{},{}\n", rna_id, pos_1 - 1, pos_2 - 1)
      });
    }
    buf_4_writer_2_output_file.push_str(&buf_4_col_pair);
  } */
  /* println!("Computed MEA CSS.");
  for pos_pairs in &mea_css.pos_pair_seqs {
    if pos_pairs[0].0 == 0. {continue;}
    let mut buf_4_col_pair = String::new();
    for (i, &(j, k)) in pos_pairs.iter().enumerate() {
      if is_gap_pos(j) || is_gap_pos(k) {continue;}
      let rna_id = mea_css.rna_ids[i];
      buf_4_col_pair.push_str(&if i < num_of_fasta_records - 1 {
        format!("{}:{},{} ", rna_id, j as Pos - 1, k as Pos - 1)
      } else {
        format!("{}:{},{}\n", rna_id, j as Pos - 1, k as Pos - 1)
      });
    }
    buf_4_writer_2_output_file.push_str(&buf_4_col_pair);
  }
  let _ = writer_2_output_file.write_all(buf_4_writer_2_output_file.as_bytes()); */
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
  let reader_2_clustal_file = BufReader::new(File::open(clustal_file_path).expect("Failed to read a CLUSTAL file."));
  // let mut buf_4_reader_2_clustal_file = Vec::new();
  let mut seq_pointer = 0;
  let mut pos_pointer = 0;
  let mut are_seq_ids_read = false;
  for (i, string) in reader_2_clustal_file.lines().enumerate() {
    let string = string.expect("Failed to read a CLUSTAL file.");
    // println!("{}", sa.cols.len());
    if i == 0 || string.len() == 0 || string.starts_with(" ") {
      // println!("len 0.");
      if sa.cols.len() > 0 {
        seq_pointer = 0;
        pos_pointer = sa.cols.len();
        are_seq_ids_read = true;
      }
      continue;
    }
    let mut substrings = string.split_whitespace();
    let substring = substrings.next().expect("Failed to read a CLUSTAL file.");
    if !are_seq_ids_read {
      seq_ids.push(String::from(substring));
    }
    let substring = substrings.next().expect("Failed to read a CLUSTAL file.");
    // println!("{}.", &substring);
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
  let num_of_rnas = sa.cols[0].len();
  let pseudo_bases = vec![PSEUDO_BASE; num_of_rnas];
  sa.cols.insert(0, pseudo_bases.clone());
  sa.cols.push(pseudo_bases);
  (sa, seq_ids)
}

#[inline]
fn get_mean_upp_mat(upp_mats_with_rna_ids: &ProbsWithRnaIds, sa: &SeqAlign) -> Probs {
  let sa_len = sa.cols.len();
  let num_of_rnas = sa.cols[0].len();
  let mut mean_upp_mat = vec![0.; sa_len];
  for i in 0 .. sa_len {
    let mut mean_upp = 0.;
    for j in 0 .. num_of_rnas {
      if sa.cols[i][j] == GAP {continue;}
      mean_upp += upp_mats_with_rna_ids[j][sa.pos_map_sets[i][j]];
    }
    mean_upp_mat[i] = mean_upp / num_of_rnas as Prob;
  }
  mean_upp_mat
}

/* #[inline]
fn get_mea_css_of_node(guide_tree: &GuideTree, node: &NodeIndex<usize>, mea_csss: &MeaCsss, gamma_plus_1: Prob, bpap_mats_with_rna_id_pairs: &Prob4dMatsWithRnaIdPairs, upp_mats_with_rna_ids: &ProbSeqsWithRnaIds, num_of_fasta_records: usize, num_of_iterative_refinements: usize) -> MeaCss {
  let node_index = node.index();
  if node.index() < num_of_fasta_records {
    mea_csss[node_index].clone()
  } else {
    let mut child_node_pair = guide_tree.neighbors_directed(*node, Outgoing);
    let child_node_pair = (
      child_node_pair.next().expect("Failed to get a child node of a node in a guide tree."),
      child_node_pair.next().expect("Failed to get a child node of a node in a guide tree."),
    );
    let mut mea_css = get_mea_consensus_ss(
      &(
        &get_mea_css_of_node(guide_tree, &child_node_pair.0, mea_csss, gamma_plus_1, bpap_mats_with_rna_id_pairs, upp_mats_with_rna_ids, num_of_fasta_records, num_of_iterative_refinements),
        &get_mea_css_of_node(guide_tree, &child_node_pair.1, mea_csss, gamma_plus_1, bpap_mats_with_rna_id_pairs, upp_mats_with_rna_ids, num_of_fasta_records, num_of_iterative_refinements),
      ),
      gamma_plus_1,
      bpap_mats_with_rna_id_pairs,
      upp_mats_with_rna_ids,
    );
    if child_node_pair.0.index() >= num_of_fasta_records || child_node_pair.1.index() >= num_of_fasta_records {
      let num_of_rnas = mea_css.rna_ids.len();
      let mut rand_num_generator = rand::thread_rng();
      for _ in 0 .. num_of_iterative_refinements {
        let mut new_mea_css = mea_css.clone();
        let index_4_remove = rand_num_generator.gen_range(0, num_of_rnas);
        let chosen_rna_id = new_mea_css.rna_ids[index_4_remove].clone();
        new_mea_css.rna_ids.remove(index_4_remove);
        new_mea_css.seq_num -= 1;
        for pos_pairs in &mut new_mea_css.pos_pair_seqs {
          pos_pairs.remove(index_4_remove);
        }
        for poss in &mut new_mea_css.pos_seqs_with_cols {
          poss.remove(index_4_remove);
        }
        new_mea_css = get_mea_consensus_ss(
          &(
            &new_mea_css,
            &get_mea_css_of_node(guide_tree, &NodeIndex::new(chosen_rna_id), mea_csss, gamma_plus_1, bpap_mats_with_rna_id_pairs, upp_mats_with_rna_ids, num_of_fasta_records, num_of_iterative_refinements),
          ),
          gamma_plus_1,
          bpap_mats_with_rna_id_pairs,
          upp_mats_with_rna_ids,
        );
        if mea_css.mea < new_mea_css.mea {
          mea_css = new_mea_css;
          println!("Refined MEA CSS.");
        }
      }
    }
    mea_css
  }
} */

#[inline]
fn get_mea_css_str(mea_css: &MeaCss, sa_len: usize) -> MeaCssStr {
  let mut mea_css_str = vec![UNPAIRING_BASE; sa_len - 2];
  let pseudo_pos_pair = (0, sa_len - 1);
  for bpa_pos_pair in &mea_css.bpa_pos_pairs {
    let (i, j) = bpa_pos_pair;
    if *bpa_pos_pair != pseudo_pos_pair {
      mea_css_str[i - 1] = BASE_PAIRING_LEFT_BASE;
      mea_css_str[j - 1] = BASE_PAIRING_RIGHT_BASE;
    }
  }
  mea_css_str
}
