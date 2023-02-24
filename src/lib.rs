extern crate consprob;

pub use consprob::*;
pub use std::fs::remove_file;
use std::process::{Command, Output};

pub type SparseScores<T> = HashMap<T, Score>;
pub type ScoresHashedPoss<T> = HashMap<T, SparseScores<T>>;
pub type ColSetsHashedCols<T> = HashMap<T, SparseScores<T>>;
pub type ColsHashedCols<T> = HashMap<T, T>;
pub type SparsePosMat<T> = HashSet<PosPair<T>>;
pub type SparseProbMats<T> = Vec<SparseProbMat<T>>;

pub const GAP: Char = b'-';
pub const EXAMPLE_CLUSTAL_FILE_PATH: &str = "assets/test_seqs.aln";
pub const MIN_LOG_HYPERPARAM: i32 = -4;
pub const MAX_LOG_HYPERPARAM: i32 = 10;

pub fn consalifold<T>(
  basepair_probs_mix: &SparseProbMat<T>,
  align: &Align<T>,
  hyperparam: Score,
) -> SparsePosMat<T>
where
  T: HashIndex,
{
  let align_len = align.cols.len();
  let align_len = T::from_usize(align_len).unwrap();
  let basepair_scores: SparseScoreMat<T> = basepair_probs_mix
    .iter()
    .filter(|x| hyperparam * x.1 - 1. >= 0.)
    .map(|x| (*x.0, *x.1))
    .collect();
  let mut scores_hashed_cols = ScoresHashedPoss::default();
  let mut right_basepairs_hashed_cols = ColSetsHashedCols::<T>::default();
  for (x, &y) in &basepair_scores {
    match right_basepairs_hashed_cols.get_mut(&x.0) {
      Some(z) => {
        z.insert(x.1, y);
      }
      None => {
        let mut z = SparseScores::default();
        z.insert(x.1, y);
        right_basepairs_hashed_cols.insert(x.0, z);
      }
    }
  }
  let mut rightmost_basepairs_hashed_cols = ColsHashedCols::<T>::default();
  for (&x, y) in &right_basepairs_hashed_cols {
    let y = y.keys().copied().max().unwrap();
    rightmost_basepairs_hashed_cols.insert(x, y);
  }
  for i in range_inclusive(T::one(), align_len).rev() {
    if let Some(&j) = rightmost_basepairs_hashed_cols.get(&i) {
      let j = (i, j);
      let x = get_scores(&scores_hashed_cols, &j);
      update_scores_hashed_cols(
        &mut scores_hashed_cols,
        j.0,
        &x,
        &right_basepairs_hashed_cols,
      );
    }
  }
  let pseudo_col_pair = (T::zero(), align_len + T::one());
  let mut basepairs = SparsePosMat::<T>::default();
  traceback(&mut basepairs, &pseudo_col_pair, &scores_hashed_cols);
  basepairs
}

pub fn traceback<T>(
  x: &mut SparsePosMat<T>,
  y: &PosPair<T>,
  z: &ScoresHashedPoss<T>,
) where
  T: HashIndex,
{
  let mut a;
  let b = get_scores(z, y);
  let (i, j) = *y;
  let mut k = j - T::one();
  while k > i {
    a = b[&k];
    let c = b[&(k - T::one())];
    if c == a {
      k = k - T::one();
    }
    if let Some(c) = z.get(&k) {
      for (&c, d) in c {
        if i >= c {
          continue;
        }
        let e = c - T::one();
        let b = b[&e];
        let d = b + d;
        if d == a {
          let d = (c, k);
          traceback(x, &d, z);
          let d = (c - T::one(), k - T::one());
          x.insert(d);
          k = e;
          break;
        }
      }
    }
  }
}

pub fn update_scores_hashed_cols<T>(
  x: &mut ScoresHashedPoss<T>,
  y: T,
  z: &SparseScores<T>,
  a: &ColSetsHashedCols<T>,
) where
  T: HashIndex,
{
  let a = &a[&y];
  for (&a, &b) in a {
    let b = b + z[&(a - T::one())];
    match x.get_mut(&a) {
      Some(x) => {
        x.insert(y, b);
      }
      None => {
        let mut z = SparseScores::default();
        z.insert(y, b);
        x.insert(a, z);
      }
    }
  }
}

pub fn get_scores<T>(x: &ScoresHashedPoss<T>, y: &PosPair<T>) -> SparseScores<T>
where
  T: HashIndex,
{
  let (i, j) = *y;
  let mut z = SparseScores::<T>::default();
  for k in range(i, j) {
    if k == i {
      z.insert(k, 0.);
      continue;
    }
    let mut y = z[&(k - T::one())];
    if let Some(x) = x.get(&k) {
      for (&l, x) in x {
        if i >= l {
          continue;
        }
        let z = z[&(l - T::one())];
        let z = x + z;
        if z > y {
          y = z;
        }
      }
    }
    z.insert(k, y);
  }
  z
}

pub fn get_basepair_probs_alifold<T>(align_file_path: &Path, output_dir_path: &Path) -> SparseProbMat<T>
where
  T: HashIndex,
{
  let cwd = env::current_dir().unwrap();
  let align_file_path = cwd.join(align_file_path);
  let align_file_prefix = align_file_path.file_stem().unwrap().to_str().unwrap();
  let arg = format!("--id-prefix={align_file_prefix}");
  let args = vec![
    "-p",
    align_file_path.to_str().unwrap(),
    &arg,
    "--noPS",
    "--noDP",
  ];
  let _ = env::set_current_dir(output_dir_path);
  let _ = run_command("RNAalifold", &args, "Failed to run RNAalifold");
  let _ = env::set_current_dir(cwd);
  let mut basepair_probs_alifold = SparseProbMat::<T>::default();
  let output_file_path = output_dir_path.join(String::from(align_file_prefix) + "_0001_ali.out");
  let output_file = BufReader::new(File::open(output_file_path.clone()).unwrap());
  for (x, y) in output_file.lines().enumerate() {
    if x == 0 {
      continue;
    }
    let y = y.unwrap();
    if !y.starts_with(' ') {
      continue;
    }
    let y: Vec<&str> = y.split_whitespace().collect();
    let z = (
      T::from_usize(y[0].parse().unwrap()).unwrap(),
      T::from_usize(y[1].parse().unwrap()).unwrap(),
    );
    let mut a = String::from(y[3]);
    a.pop();
    let a = 0.01 * a.parse::<Prob>().unwrap();
    if a == 0. {
      continue;
    }
    basepair_probs_alifold.insert(z, a);
  }
  let _ = remove_file(output_file_path);
  basepair_probs_alifold
}

pub fn get_basepair_probs_mix<T>(
  align: &Align<T>,
  basepair_prob_mats: &SparseProbMats<T>,
  basepair_probs_alifold: &SparseProbMat<T>,
  mix_weight: Prob,
) -> SparseProbMat<T>
where
  T: HashIndex,
{
  let mut basepair_probs_mix = SparseProbMat::<T>::default();
  let align_len = align.cols.len();
  let num_rnas = align.cols[0].len();
  for i in 0..align_len {
    let pos_maps = &align.pos_map_sets[i];
    let short_i = T::from_usize(i).unwrap();
    for j in i + 1..align_len {
      let short_j = T::from_usize(j).unwrap();
      let pos_pair = (short_i, short_j);
      let basepair_prob_alifold = match basepair_probs_alifold.get(&pos_pair) {
        Some(&x) => x,
        None => 0.,
      };
      let pos_maps2 = &align.pos_map_sets[j];
      let pos_map_pairs: Vec<(T, T)> = pos_maps
        .iter()
        .zip(pos_maps2.iter())
        .map(|(&x, &y)| (x, y))
        .collect();
      let mut basepair_prob_sum = 0.;
      for (x, y) in pos_map_pairs.iter().zip(basepair_prob_mats.iter()) {
        if let Some(&y) = y.get(x) {
          basepair_prob_sum += y;
        }
      }
      let basepair_prob_avg = basepair_prob_sum / num_rnas as Prob;
      let basepair_prob_mix = mix_weight * basepair_prob_avg + (1. - mix_weight) * basepair_prob_alifold;
      let pos_pair = (pos_pair.0 + T::one(), pos_pair.1 + T::one());
      basepair_probs_mix.insert(pos_pair, basepair_prob_mix);
    }
  }
  basepair_probs_mix
}

pub fn base2char(c: Base) -> Char {
  match c {
    A => A_UPPER,
    C => C_UPPER,
    G => G_UPPER,
    U => U_UPPER,
    PSEUDO_BASE => GAP,
    _ => {
      panic!();
    }
  }
}

pub fn run_command(x: &str, y: &[&str], z: &str) -> Output {
  Command::new(x).args(y).output().expect(z)
}
