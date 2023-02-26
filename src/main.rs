extern crate consalifold;
extern crate num_cpus;

use consalifold::*;
use std::env;
use std::fs::create_dir;
use std::fs::remove_file;
use std::fs::File;
use std::io::{BufRead, BufWriter};
use std::path::Path;

type InputsConsalifoldMultithreaded<'a> = (
  &'a mut Pool,
  &'a Cols,
  &'a SeqIds,
  Prob,
  Prob,
  &'a Path,
  Score,
  Prob,
  &'a Path,
  ScoreModel,
);

const DEFAULT_MIX_WEIGHT: Prob = 0.5;
const DEFAULT_HYPERPARAM: Prob = NEG_INFINITY;
enum ScoreModel {
  Turner,
  Posterior,
}
const DEFAULT_SCORE_MODEL: &str = "turner";
const README_CONTENTS_ALIFOLD: &str = "# hyperparam=x.sth\nThis file type contains a predicted consensus secondary structure in Stockholm format.\n\n";

fn main() {
  let args = env::args().collect::<Vec<Arg>>();
  let program_name = args[0].clone();
  let mut opts = Options::new();
  opts.reqopt(
    "i",
    "input_file_path",
    "An input CLUSTAL/FASTA/STOCKHOLM file path containing an RNA alignment",
    "STR",
  );
  opts.reqopt("o", "output_dir_path", "An output directory path", "STR");
  opts.optopt(
    "",
    "min_basepair_prob",
    &format!("A minimum base-pairing probability (Use {DEFAULT_MIN_BASEPAIR_PROB} by default)"),
    "FLOAT",
  );
  opts.optopt(
    "",
    "min_match_prob",
    &format!("A minimum matching probability (Use {DEFAULT_MIN_MATCH_PROB} by default)"),
    "FLOAT",
  );
  opts.optopt(
    "p",
    "hyperparam",
    "A specific hyper-parameter rather than a range of hyper-parameters",
    "FLOAT",
  );
  opts.optopt(
    "",
    "mix_weight",
    &format!("A mixture weight (Use {DEFAULT_MIX_WEIGHT} by default)"),
    "FLOAT",
  );
  opts.optopt(
    "m",
    "score_model",
    &format!(
      "Choose a structural alignment scoring model from turner, posterior (Use {DEFAULT_SCORE_MODEL} by default)"
    ),
    "STR",
  );
  opts.optopt(
    "t",
    "num_threads",
    "The number of threads in multithreading (Use all the threads of this computer by default)",
    "UINT",
  );
  opts.optflag("h", "help", "Print a help menu");
  let matches = match opts.parse(&args[1..]) {
    Ok(opt) => opt,
    Err(failure) => {
      print_program_usage(&program_name, &opts);
      panic!("{}", failure.to_string())
    }
  };
  if matches.opt_present("h") {
    print_program_usage(&program_name, &opts);
    return;
  }
  let input_file_path = matches.opt_str("i").unwrap();
  let input_file_path = Path::new(&input_file_path);
  let output_dir_path = matches.opt_str("o").unwrap();
  let output_dir_path = Path::new(&output_dir_path);
  let score_model = if matches.opt_present("m") {
    let x = matches.opt_str("m").unwrap();
    if x == "turner" {
      ScoreModel::Turner
    } else if x == "posterior" {
      ScoreModel::Posterior
    } else {
      panic!();
    }
  } else {
    ScoreModel::Turner
  };
  let min_basepair_prob = if matches.opt_present("min_basepair_prob") {
    matches
      .opt_str("min_basepair_prob")
      .unwrap()
      .parse()
      .unwrap()
  } else {
    DEFAULT_MIN_BASEPAIR_PROB
  };
  let min_match_prob = if matches.opt_present("min_match_prob") {
    matches.opt_str("min_match_prob").unwrap().parse().unwrap()
  } else {
    DEFAULT_MIN_MATCH_PROB
  };
  let hyperparam = if matches.opt_present("hyperparam") {
    matches.opt_str("hyperparam").unwrap().parse().unwrap()
  } else {
    DEFAULT_HYPERPARAM
  };
  let mix_weight = if matches.opt_present("mix_weight") {
    matches.opt_str("mix_weight").unwrap().parse().unwrap()
  } else {
    DEFAULT_MIX_WEIGHT
  };
  let num_threads = if matches.opt_present("t") {
    matches.opt_str("t").unwrap().parse().unwrap()
  } else {
    num_cpus::get() as NumThreads
  };
  let extension = input_file_path.extension().unwrap().to_str().unwrap();
  let matches_stockholm = matches!(extension, "sto" | "stk" | "sth");
  let matches_fasta = matches!(extension, "fasta" | "fna" | "ffn" | "faa" | "frn" | "fa");
  let (cols, seq_ids) = if matches_stockholm {
    read_align_stockholm(input_file_path)
  } else if matches_fasta {
    read_align_fasta(input_file_path)
  } else {
    read_align_clustal(input_file_path)
  };
  let align_len = cols.len();
  let mut thread_pool = Pool::new(num_threads);
  if align_len + 2 <= u8::MAX as usize {
    consalifold_multithreaded::<u8>((
      &mut thread_pool,
      &cols,
      &seq_ids,
      min_basepair_prob,
      min_match_prob,
      output_dir_path,
      hyperparam,
      mix_weight,
      input_file_path,
      score_model,
    ));
  } else {
    consalifold_multithreaded::<u16>((
      &mut thread_pool,
      &cols,
      &seq_ids,
      min_basepair_prob,
      min_match_prob,
      output_dir_path,
      hyperparam,
      mix_weight,
      input_file_path,
      score_model,
    ));
  }
}

fn consalifold_multithreaded<T>(inputs: InputsConsalifoldMultithreaded)
where
  T: HashIndex,
{
  let (
    thread_pool,
    cols,
    seq_ids,
    min_basepair_prob,
    min_match_prob,
    output_dir_path,
    hyperparam,
    mix_weight,
    input_file_path,
    score_model,
  ) = inputs;
  let matches_posterior_model = matches!(score_model, ScoreModel::Posterior);
  let mut align = Align::<T>::new();
  align.cols = cols.clone();
  let num_rnas = align.cols[0].len();
  let mut seq_lens = vec![0_usize; num_rnas];
  let align_len = align.cols.len();
  align.pos_map_sets = vec![vec![T::zero(); num_rnas]; align_len];
  let mut fasta_records = vec![FastaRecord::origin(); num_rnas];
  for i in 0..align_len {
    for j in 0..num_rnas {
      let x = align.cols[i][j];
      if x != PSEUDO_BASE {
        fasta_records[j].seq.push(x);
        seq_lens[j] += 1;
        align.pos_map_sets[i][j] = T::from_usize(seq_lens[j]).unwrap();
      }
    }
  }
  for i in 0..num_rnas {
    fasta_records[i].seq.insert(0, PSEUDO_BASE);
    fasta_records[i].seq.push(PSEUDO_BASE);
    fasta_records[i].fasta_id = seq_ids[i].clone();
  }
  if !output_dir_path.exists() {
    let _ = create_dir(output_dir_path);
  }
  let mut align_scores = AlignScores::new(0.);
  align_scores.transfer();
  let seqs = fasta_records.iter().map(|x| &x.seq[..]).collect();
  let basepair_probs_alifold = get_basepair_probs_alifold(input_file_path, output_dir_path);
  let produces_context_profs = false;
  let produces_match_probs = false;
  let alignfold_prob_mats_avg = if matches_posterior_model {
    locarnap(thread_pool, &fasta_records, output_dir_path)
  } else {
    let x = consprob::<T>(
      thread_pool,
      &seqs,
      min_basepair_prob,
      min_match_prob,
      produces_context_profs,
      produces_match_probs,
      &align_scores,
    ).0;
    write_alignfold_prob_mats::<T>(
      output_dir_path,
      &x,
      &MatchProbsHashedIds::<T>::default(),
      produces_context_profs,
      produces_match_probs,
    );
    x
  };
  let basepair_prob_mats = alignfold_prob_mats_avg.iter().map(|x| x.basepair_probs.clone()).collect();
  let basepair_probs_mix = get_basepair_probs_mix(&align, &basepair_prob_mats, &basepair_probs_alifold, mix_weight);
  if hyperparam != NEG_INFINITY {
    let output_file_path = output_dir_path.join(format!("hyperparam={hyperparam}.sth"));
    let hyperparam = hyperparam + 1.;
    consalifold_wrapped(
      &basepair_probs_mix,
      &align,
      hyperparam,
      &output_file_path,
      &fasta_records,
    );
  } else {
    thread_pool.scoped(|x| {
      let y = &basepair_probs_mix;
      let z = &align;
      let a = &fasta_records;
      for b in MIN_LOG_HYPERPARAM..MAX_LOG_HYPERPARAM + 1 {
        let b = (2. as Score).powi(b);
        let output_file_path = output_dir_path.join(format!("hyperparam={b}.sth"));
        let b = b + 1.;
        x.execute(move || {
          consalifold_wrapped::<T>(y, z, b, &output_file_path, a);
        });
      }
    });
  }
  let mut readme_contents = String::from(README_CONTENTS_ALIFOLD);
  readme_contents.push_str(README_CONTENTS);
  write_readme(output_dir_path, &readme_contents);
}

fn locarnap<T>(
  thread_pool: &mut Pool,
  fasta_records: &FastaRecords,
  output_dir_path: &Path,
) -> ProbMatSetsAvg<T>
where
  T: HashIndex,
{
  let num_fasta_records = fasta_records.len();
  for (i, fasta_record) in fasta_records.iter().enumerate() {
    let output_file_path = output_dir_path.join(format!("locarnap_seq{i}.fa"));
    let mut writer = BufWriter::new(File::create(output_file_path).unwrap());
    let seq = &fasta_record.seq;
    let seq_len = seq.len();
    let buf = format!(">{}\n{}", i, seq2str(&seq[1..seq_len - 1]));
    let _ = writer.write_all(buf.as_bytes());
  }
  let mut alignfold_probs_hashed_ids = AlignfoldProbsHashedIds::<T>::default();
  for x in 0..num_fasta_records {
    for y in x + 1..num_fasta_records {
      let y = (x, y);
      alignfold_probs_hashed_ids.insert(y, AlignfoldProbMats::<T>::origin());
    }
  }
  thread_pool.scoped(|x| {
    for (y, z) in alignfold_probs_hashed_ids.iter_mut() {
      x.execute(move || {
        z.basepair_probs_pair = exec_locarnap(y, output_dir_path);
      });
    }
  });
  let mut alignfold_prob_mats_avg = vec![AlignfoldProbMatsAvg::<T>::origin(); num_fasta_records];
  thread_pool.scoped(|x| {
    let y = &alignfold_probs_hashed_ids;
    for (z, a) in alignfold_prob_mats_avg.iter_mut().enumerate() {
      let b = fasta_records[z].seq.len();
      let c = output_dir_path.join(format!("locarnap_seq{z}.fa"));
      x.execute(move || {
        *a = postprocess_locarnap::<T>(y, z, b, num_fasta_records);
        let _ = remove_file(c);
      });
    }
  });
  alignfold_prob_mats_avg
}

fn exec_locarnap<T>(rna_id_pair: &RnaIdPair, output_dir_path: &Path) -> SparseProbMatPair<T>
where
  T: HashIndex,
{
  let mut basepair_probs_pair = (SparseProbMat::<T>::default(), SparseProbMat::<T>::default());
  let (seq_file_path, seq_file_path2) = (
    output_dir_path.join(format!("locarnap_seq{}.fa", rna_id_pair.0)),
    output_dir_path.join(format!("locarnap_seq{}.fa", rna_id_pair.1)),
  );
  let output_file_path = output_dir_path.join(format!(
    "locarnap_seq{}_seq{}.dat",
    rna_id_pair.0, rna_id_pair.1
  ));
  let arg = format!(
    "--write-arcmatch-probs={}",
    output_file_path.to_str().unwrap()
  );
  let args = vec![
    seq_file_path.to_str().unwrap(),
    seq_file_path2.to_str().unwrap(),
    &arg,
  ];
  let _ = run_command("locarna_p", &args, "Failed to run LocARNA-P");
  let output_file = BufReader::new(File::open(output_file_path.clone()).unwrap());
  for x in output_file.lines() {
    let x: Vec<String> = x.unwrap().split_whitespace().map(String::from).collect();
    let (i, j, k, l, y) = (
      T::from_usize(x[0].parse().unwrap()).unwrap(),
      T::from_usize(x[1].parse().unwrap()).unwrap(),
      T::from_usize(x[2].parse().unwrap()).unwrap(),
      T::from_usize(x[3].parse().unwrap()).unwrap(),
      x[4].parse().unwrap(),
    );
    match basepair_probs_pair.0.get_mut(&(i, j)) {
      Some(x) => {
        *x += y;
      }
      None => {
        basepair_probs_pair.0.insert((i, j), y);
      }
    }
    match basepair_probs_pair.1.get_mut(&(k, l)) {
      Some(x) => {
        *x += y;
      }
      None => {
        basepair_probs_pair.1.insert((k, l), y);
      }
    }
  }
  let _ = remove_file(output_file_path);
  basepair_probs_pair
}

pub fn seq2str(x: &[Base]) -> String {
  let mut y = Vec::<u8>::new();
  for &x in x {
    let x = match x {
      A => A_UPPER,
      C => C_UPPER,
      G => G_UPPER,
      U => U_UPPER,
      _ => {
        panic!();
      }
    };
    y.push(x);
  }
  String::from_utf8(y).unwrap()
}

fn postprocess_locarnap<T>(
  alignfold_probs_hashed_ids: &AlignfoldProbsHashedIds<T>,
  rna_id: RnaId,
  unpair_probs_len: usize,
  num_rnas: usize,
) -> AlignfoldProbMatsAvg<T>
where
  T: HashIndex,
{
  let weight = 1. / (num_rnas - 1) as Prob;
  let mut alignfold_prob_mats_avg = AlignfoldProbMatsAvg::new(unpair_probs_len);
  for x in 0..num_rnas {
    if x == rna_id {
      continue;
    }
    let y = if rna_id < x { (rna_id, x) } else { (x, rna_id) };
    let y = &alignfold_probs_hashed_ids[&y];
    let y = if rna_id < x {
      &y.basepair_probs_pair.0
    } else {
      &y.basepair_probs_pair.1
    };
    for (x, &y) in y.iter() {
      let y = weight * y;
      match alignfold_prob_mats_avg.basepair_probs.get_mut(x) {
        Some(z) => {
          *z += y;
        }
        None => {
          alignfold_prob_mats_avg.basepair_probs.insert(*x, y);
        }
      }
    }
  }
  alignfold_prob_mats_avg
}

fn consalifold_wrapped<T>(
  basepair_probs_mix: &SparseProbMat<T>,
  align: &Align<T>,
  hyperparam: Prob,
  output_file_path: &Path,
  fasta_records: &FastaRecords,
) where
  T: HashIndex,
{
  let basepairs = consalifold::<T>(basepair_probs_mix, align, hyperparam);
  let mut writer = BufWriter::new(File::create(output_file_path).unwrap());
  let mut buf = "# STOCKHOLM 1.0\n".to_string();
  let align_len = align.cols.len();
  let max_seq_id_len = fasta_records
    .iter()
    .map(|x| x.fasta_id.len())
    .max()
    .unwrap();
  let descriptor = "#=GC SS_cons";
  let descriptor_len = descriptor.len();
  let max_seq_id_len = max_seq_id_len.max(descriptor_len);
  for (x, y) in fasta_records.iter().enumerate() {
    let y = &y.fasta_id;
    buf.push_str(y);
    let mut y = vec![b' '; max_seq_id_len - y.len() + 2];
    let mut z = (0..align_len)
      .map(|y| base2char(align.cols[y][x]))
      .collect::<Vec<Char>>();
    y.append(&mut z);
    let y = unsafe { from_utf8_unchecked(&y) };
    buf.push_str(y);
    buf.push('\n');
  }
  buf.push_str(descriptor);
  let mut stockholm_row = vec![b' '; max_seq_id_len - descriptor_len + 2];
  let mut fold_str = get_fold_str(&basepairs, align_len);
  stockholm_row.append(&mut fold_str);
  let stockholm_row = unsafe { from_utf8_unchecked(&stockholm_row) };
  buf.push_str(stockholm_row);
  buf.push_str("\n//");
  let _ = writer.write_all(buf.as_bytes());
}

fn get_fold_str<T>(x: &SparsePosMat<T>, y: usize) -> FoldStr
where
  T: HashIndex,
{
  let mut z = vec![UNPAIR; y];
  for &(i, j) in x {
    z[i.to_usize().unwrap()] = BASEPAIR_LEFT;
    z[j.to_usize().unwrap()] = BASEPAIR_RIGHT;
  }
  z
}
