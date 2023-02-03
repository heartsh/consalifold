extern crate consalifold;
extern crate criterion;

use consalifold::*;
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_consalifold(criterion: &mut Criterion) {
  let num_of_threads = num_cpus::get() as NumOfThreads;
  let mut thread_pool = Pool::new(num_of_threads);
  let input_file_path = Path::new(EXAMPLE_CLUSTAL_FILE_PATH);
  let (cols, seq_ids) = read_sa_from_clustal_file(input_file_path);
  let mut sa = SeqAlign::<u8>::new();
  sa.cols = cols.clone();
  let num_of_rnas = sa.cols[0].len();
  let mut seq_lens = vec![0 as usize; num_of_rnas];
  let sa_len = sa.cols.len();
  sa.pos_map_sets = vec![vec![0; num_of_rnas]; sa_len];
  let mut fasta_records = vec![FastaRecord::origin(); num_of_rnas];
  for i in 0..sa_len {
    for j in 0..num_of_rnas {
      let base = sa.cols[i][j];
      if base != PSEUDO_BASE {
        fasta_records[j].seq.push(base);
        seq_lens[j] += 1;
        sa.pos_map_sets[i][j] = u8::from_usize(seq_lens[j]).unwrap();
      }
    }
  }
  for i in 0..num_of_rnas {
    fasta_records[i].seq.insert(0, PSEUDO_BASE);
    fasta_records[i].seq.push(PSEUDO_BASE);
    fasta_records[i].fasta_id = seq_ids[i].clone();
  }
  let mut align_feature_score_sets = AlignFeatureCountSets::new(0.);
  align_feature_score_sets.transfer();
  let seqs = fasta_records.iter().map(|x| &x.seq[..]).collect();
  let prob_mat_sets = consprob::<u8>(
      &mut thread_pool,
      &seqs,
      DEFAULT_MIN_BPP,
      DEFAULT_MIN_ALIGN_PROB,
      false,
      false,
      &align_feature_score_sets,
    )
    .0;
  let bpp_mats = prob_mat_sets.iter().map(|x| x.bpp_mat.clone()).collect();
  let rnaalifold_bpp_mat = SparseProbMat::<u8>::default();
  let mix_bpp_mat = get_mix_bpp_mat(&sa, &bpp_mats, &rnaalifold_bpp_mat, 1.);
  let gamma = (2. as Prob).powi(MAX_POW_OF_2) + 1.;
  criterion.bench_function("consalifold::<u8>", |b| {
    b.iter(|| {
      let _ = consalifold::<u8>(&mix_bpp_mat, &sa, gamma);
    });
  });
}

criterion_group!(benches, bench_consalifold);
criterion_main!(benches);
