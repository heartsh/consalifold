extern crate consalifold;
extern crate criterion;

use consalifold::*;
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_consalifold(criterion: &mut Criterion) {
  let num_threads = num_cpus::get() as NumThreads;
  let mut thread_pool = Pool::new(num_threads);
  let input_file_path = Path::new(EXAMPLE_CLUSTAL_FILE_PATH);
  let (cols, seq_ids) = read_align_clustal(input_file_path);
  let mut align = Align::<u8>::new();
  align.cols = cols;
  let num_rnas = align.cols[0].len();
  let mut seq_lens = vec![0_usize; num_rnas];
  let align_len = align.cols.len();
  align.pos_map_sets = vec![vec![0; num_rnas]; align_len];
  let mut fasta_records = vec![FastaRecord::origin(); num_rnas];
  for i in 0..align_len {
    for j in 0..num_rnas {
      let x = align.cols[i][j];
      if x != PSEUDO_BASE {
        fasta_records[j].seq.push(x);
        seq_lens[j] += 1;
        align.pos_map_sets[i][j] = u8::from_usize(seq_lens[j]).unwrap();
      }
    }
  }
  for i in 0..num_rnas {
    fasta_records[i].seq.insert(0, PSEUDO_BASE);
    fasta_records[i].seq.push(PSEUDO_BASE);
    fasta_records[i].fasta_id = seq_ids[i].clone();
  }
  let mut align_scores = AlignScores::new(0.);
  align_scores.transfer();
  let seqs = fasta_records.iter().map(|x| &x.seq[..]).collect();
  let produces_struct_profs = false;
  let produces_match_probs = false;
  let alignfold_prob_mats_avg = consprob::<u8>(
    &mut thread_pool,
    &seqs,
    DEFAULT_MIN_BASEPAIR_PROB,
    DEFAULT_MIN_MATCH_PROB,
    produces_struct_profs,
    produces_match_probs,
    &align_scores,
  )
  .0;
  let basepair_prob_mats = alignfold_prob_mats_avg
    .iter()
    .map(|x| x.basepair_probs.clone())
    .collect();
  let basepair_probs_alifold = SparseProbMat::<u8>::default();
  let mix_weight = 1.;
  let basepair_probs_mix = get_basepair_probs_mix(
    &align,
    &basepair_prob_mats,
    &basepair_probs_alifold,
    mix_weight,
  );
  let hyperparam = (2. as Score).powi(MAX_LOG_HYPERPARAM) + 1.;
  criterion.bench_function("consalifold::<u8>", |x| {
    x.iter(|| {
      let _ = consalifold::<u8>(&basepair_probs_mix, &align, hyperparam);
    });
  });
}

criterion_group!(benches, bench_consalifold);
criterion_main!(benches);
