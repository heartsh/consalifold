extern crate mearcof;
extern crate getopts;
extern crate scoped_threadpool;
extern crate bio;
extern crate num_cpus;
extern crate petgraph;
extern crate rand;

use mearcof::*;
use getopts::Options;
use self::scoped_threadpool::Pool;
use std::env;
use std::path::Path;
use bio::io::fasta::Reader;
use std::io::prelude::*;
use std::io::{BufReader, BufWriter};
use std::fs::File;
use petgraph::{Outgoing, Incoming};
use rand::Rng;

type Arg = String;
type NumOfThreads = u32;
type FastaId = String;
type FastaRecord = (FastaId, Seq, usize);
type FastaRecords = Vec<FastaRecord>;
type Strings = Vec<String>;
type MeaCsss = Vec<MeaCss>;

const DEFAULT_GAMMA: Prob = 9.;
const DEFAULT_NUM_OF_ITERATIVE_REFINEMENTS: usize = 5;
const VERSION: &'static str = "0.1.0";

fn main() {
  let args = env::args().collect::<Vec<Arg>>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt("f", "input_fasta_file_path", "The path to an input FASTA file containing RNA sequences", "STR");
  opts.reqopt("p", "input_base_pair_align_prob_matrix_file_path", "The path to an input file containing base pair alignment probability matrices", "STR");
  opts.reqopt("o", "output_file_path", "The path to an output file which will contain estimated consensus secondary structures", "STR");
  opts.optopt("", "gamma", &format!("An MEA gamma (Uses {} by default)", DEFAULT_GAMMA), "FLOAT");
  opts.optopt("", "num_of_iterative_refinements", &format!("The number of iterative refinements (Uses {} by default)", DEFAULT_NUM_OF_ITERATIVE_REFINEMENTS), "UINT");
  opts.optopt("t", "num_of_threads", "The number of threads in multithreading (Uses the number of all the threads of this computer by default)", "UINT");
  opts.optflag("h", "help", "Print a help menu");
  let opts = match opts.parse(&args[1 ..]) {
    Ok(opt) => {opt}
    Err(failure) => {print_program_usage(&program_name, &opts); panic!(failure.to_string())}
  };
  let input_fasta_file_path = opts.opt_str("f").expect("Failed to get the path to an input FASTA file containing RNA sequences from command arguments.");
  let input_fasta_file_path = Path::new(&input_fasta_file_path);
  let input_bpap_mat_file_path = opts.opt_str("p").expect("Failed to get the path to an input file containing base pair alignment probability matrices from command arguments.");
  let input_bpap_mat_file_path = Path::new(&input_bpap_mat_file_path);
  let output_file_path = opts.opt_str("o").expect("Failed to get the path to an output file which will contain estimated consensus secondary structures from command arguments.");
  let output_file_path = Path::new(&output_file_path);
  let gamma_plus_1 = if opts.opt_present("gamma") {
    opts.opt_str("gamma").expect("Failed to get an MEA gamma from command arguments.").parse().expect("Failed to parse an MEA gamma.")
  } else {
    DEFAULT_GAMMA
  } + 1.;
  let num_of_iterative_refinements = if opts.opt_present("num_of_iterative_refinements") {
    opts.opt_str("num_of_iterative_refinements").expect("Failed to get the number of iterative refinements from command arguments.").parse().expect("Failed to parse the number of iterative refinements.")
  } else {
    DEFAULT_NUM_OF_ITERATIVE_REFINEMENTS
  };
  let num_of_threads = if opts.opt_present("t") {
    opts.opt_str("t").expect("Failed to get the number of threads in multithreading from command arguments.").parse().expect("Failed to parse the number of threads in multithreading.")
  } else {
    num_cpus::get() as NumOfThreads
  };
  let fasta_file_reader = Reader::from_file(Path::new(&input_fasta_file_path)).expect("Failed to set a FASTA file reader.");
  let mut fasta_records = FastaRecords::new();
  for fasta_record in fasta_file_reader.records() {
    let fasta_record = fasta_record.expect("Failed to read a FASTA record.");
    let seq = unsafe {from_utf8_unchecked(fasta_record.seq()).to_uppercase().as_bytes().iter().filter(|&&base| {is_rna_base(base)}).map(|&base| {base}).collect::<Seq>()};
    let seq_len = seq.len();
    fasta_records.push((String::from(fasta_record.id()), seq, seq_len));
  }
  let num_of_fasta_records = fasta_records.len();
  let mut bpap_mats_with_rna_id_pairs = Prob4dMatsWithRnaIdPairs::default();
  for rna_id_1 in 0 .. num_of_fasta_records {
    for rna_id_2 in rna_id_1 + 1 .. num_of_fasta_records {
      bpap_mats_with_rna_id_pairs.insert((rna_id_1, rna_id_2), Prob4dMat::default());
    }
  }
  let mut reader_2_input_bpap_mat_file = BufReader::new(File::open(input_bpap_mat_file_path).expect("Failed to read an input file."));
  let mut buf_4_reader_2_input_bpap_mat_file = Vec::new();
  for _ in 0 .. 2 {
    let _ = reader_2_input_bpap_mat_file.read_until(b'>', &mut buf_4_reader_2_input_bpap_mat_file);
  }
  for (i, vec) in reader_2_input_bpap_mat_file.split(b'>').enumerate() {
    if i == num_of_fasta_records * (num_of_fasta_records - 1) / 2 {continue;}
    let vec = vec.expect("Failed to read an input file.");
    let substrings = unsafe {String::from_utf8_unchecked(vec).split_whitespace().map(|string| {String::from(string)}).collect::<Strings>()};
    let rna_id_pair = substrings[0].split(',').map(|string| {String::from(string)}).collect::<Strings>();
    let rna_id_pair = (
      rna_id_pair[0].parse::<RnaId>().expect("Failed to parse an RNA ID."),
      rna_id_pair[1].parse::<RnaId>().expect("Failed to parse an RNA ID."),
    );
    let seq_len_pair = (fasta_records[rna_id_pair.0].2, fasta_records[rna_id_pair.1].2);
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
  let mut mea_csss = vec![MeaCss::new(); num_of_fasta_records];
  for i in 0 .. num_of_fasta_records {
    let ref mut mea_css = mea_csss[i];
    mea_css.seq_num = 1;
    let seq_len = fasta_records[i].2;
    mea_css.col_num = seq_len + 2;
    for j in 0 .. seq_len + 2 {
      for k in j + 1 .. seq_len + 2 {
        let pos_pair = (j, k);
        mea_css.rna_id_pos_triple_seqs_with_col_pairs.insert(pos_pair, vec![(i, pos_pair.0, pos_pair.1)]);
      }
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
      scope.execute(move || {
        *mea_mat = get_mea_consensus_ss(&mea_css_pair, gamma_plus_1, ref_2_bpap_mats_with_rna_id_pairs).mea;
      });
    }
  });
  let guide_tree = get_guide_tree(&mea_mat, num_of_fasta_records);
  let root_node = guide_tree.externals(Incoming).next().expect("Failed to get the root node of a guide tree.");
  let mea_css = get_mea_css_of_node(&guide_tree, &root_node, &mea_csss, gamma_plus_1, &bpap_mats_with_rna_id_pairs, num_of_fasta_records, num_of_iterative_refinements);
  let mut writer_2_output_file = BufWriter::new(File::create(output_file_path).expect("Failed to create an output file."));
  let mut buf_4_writer_2_output_file = format!("; The version {} of the MEARCOF program.\n; The path to the input FASTA file for computing the consensus secondary structures (= CSSs) in this file = \"{}\".\n; The path to the input base pair alignment matrix file for computing these structures = \"{}\".\n; The values of the parameters used for computing these structures are as follows.\n; \"gamma\" = {}, \"num_of_iterative_refinements\" = {}, \"num_of_threads\" = {}.\n; Each row is with pairs of the ID of each of RNA sequences and corresponding position pairs.\n\nmea = {}\n", VERSION, input_fasta_file_path.display(), input_bpap_mat_file_path.display(), gamma_plus_1 - 1., num_of_iterative_refinements, num_of_threads, mea_css.mea);
  // let seq_num = mea_css.seq_num;
  for (col_pair, rna_id_pos_triples) in mea_css.rna_id_pos_triple_seqs_with_col_pairs.iter() {
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
  }
  let _ = writer_2_output_file.write_all(buf_4_writer_2_output_file.as_bytes());
}

#[inline]
fn print_program_usage(program_name: &str, opts: &Options) {
  let program_usage = format!("The usage of this program: {} [options]", program_name);
  print!("{}", opts.usage(&program_usage));
}

#[inline]
fn get_mea_css_of_node(guide_tree: &GuideTree, node: &NodeIndex<usize>, mea_csss: &MeaCsss, gamma_plus_1: Prob, bpap_mats_with_rna_id_pairs: &Prob4dMatsWithRnaIdPairs, num_of_fasta_records: usize, num_of_iterative_refinements: usize) -> MeaCss {
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
        &get_mea_css_of_node(guide_tree, &child_node_pair.0, mea_csss, gamma_plus_1, bpap_mats_with_rna_id_pairs, num_of_fasta_records, num_of_iterative_refinements),
        &get_mea_css_of_node(guide_tree, &child_node_pair.1, mea_csss, gamma_plus_1, bpap_mats_with_rna_id_pairs, num_of_fasta_records, num_of_iterative_refinements),
      ),
      gamma_plus_1,
      bpap_mats_with_rna_id_pairs,
    );
    if child_node_pair.0.index() >= num_of_fasta_records && child_node_pair.1.index() >= num_of_fasta_records {
      let num_of_rnas = mea_css.rna_ids.len();
      let mut rand_num_generator = rand::thread_rng();
      for _ in 0 .. num_of_iterative_refinements {
        let mut new_mea_css = mea_css.clone();
        let index_4_remove = rand_num_generator.gen_range(0, num_of_rnas);
        let chosen_rna_id = new_mea_css.rna_ids[index_4_remove].clone();
        new_mea_css.rna_ids.remove(index_4_remove);
        new_mea_css.seq_num -= 1;
        for rna_id_pos_pairs in &mut new_mea_css.rna_id_pos_pair_seqs_with_cols {
          match rna_id_pos_pairs.iter().enumerate().find(|(_, &(rna_id, _))| {rna_id == chosen_rna_id}) {
            Some((index_4_remove, _)) => {
              rna_id_pos_pairs.remove(index_4_remove);
            },
            None => {},
          }
        }
        for rna_id_pos_triples in new_mea_css.rna_id_pos_triple_seqs_with_col_pairs.values_mut() {
          match rna_id_pos_triples.iter().enumerate().find(|(_, &(rna_id, _, _))| {rna_id == chosen_rna_id}) {
            Some((index_4_remove, _)) => {
              rna_id_pos_triples.remove(index_4_remove);
            },
            None => {},
          }
        }
        new_mea_css = get_mea_consensus_ss(
          &(
            &new_mea_css,
            &get_mea_css_of_node(guide_tree, &NodeIndex::new(chosen_rna_id), mea_csss, gamma_plus_1, bpap_mats_with_rna_id_pairs, num_of_fasta_records, num_of_iterative_refinements),
          ),
          gamma_plus_1,
          bpap_mats_with_rna_id_pairs,
        );
        if mea_css.mea < new_mea_css.mea {
          mea_css = new_mea_css;
        }
      }
    }
    mea_css
  }
}
