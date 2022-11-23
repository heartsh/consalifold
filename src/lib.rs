extern crate consprob;

pub use consprob::*;
pub use std::fs::remove_file;
use std::process::{Command, Output};

pub type Prob4dMatsWithRnaIdPairs<T> = HashMap<RnaIdPair, Prob4dMat<T>>;
pub type ProbsWithRnaIds = Vec<Probs>;
pub type ProbsWithPosPairs<T> = HashMap<PosPair<T>, Prob>;
pub type ProbMatsWithRnaIds<T> = Vec<SparseProbMat<T>>;
pub type SparseProbs<T> = HashMap<T, Prob>;
pub type MeaSetsWithPoss<T> = HashMap<T, SparseProbs<T>>;
pub type ColSetsWithCols<T> = HashMap<T, SparseProbs<T>>;
pub type ColsWithCols<T> = HashMap<T, T>;
pub type SparsePosMat<T> = HashSet<PosPair<T>>;
pub type SparseProbMats<T> = Vec<SparseProbMat<T>>;

pub const GAP: Char = '-' as Char;

pub fn consalifold<T>(
  mix_bpp_mat: &SparseProbMat<T>,
  sa: &SeqAlign<T>,
  gamma: Prob,
) -> SparsePosMat<T>
where
  T: HashIndex,
{
  let sa_len = sa.cols.len();
  let sa_len = T::from_usize(sa_len).unwrap();
  let mix_bpp_mat: SparseProbMat<T> = mix_bpp_mat
    .iter()
    .filter(|x| gamma * x.1 - 1. >= 0.)
    .map(|x| (*x.0, *x.1))
    .collect();
  let mut mea_sets_with_cols = MeaSetsWithPoss::default();
  let mut right_bp_cols_with_cols = ColSetsWithCols::<T>::default();
  for (col_pair, &mix_bpp) in &mix_bpp_mat {
    match right_bp_cols_with_cols.get_mut(&col_pair.0) {
      Some(cols) => {
        cols.insert(col_pair.1, mix_bpp);
      }
      None => {
        let mut cols = SparseProbs::default();
        cols.insert(col_pair.1, mix_bpp);
        right_bp_cols_with_cols.insert(col_pair.0, cols);
      }
    }
  }
  let mut rightmost_bp_cols_with_cols = ColsWithCols::<T>::default();
  for (&i, cols) in &right_bp_cols_with_cols {
    let max = cols.keys().map(|&x| x).max().unwrap();
    rightmost_bp_cols_with_cols.insert(i, max);
  }
  for i in range_inclusive(T::one(), sa_len).rev() {
    match rightmost_bp_cols_with_cols.get(&i) {
      Some(&j) => {
        let col_pair = (i, j);
        let meas = get_meas(&mea_sets_with_cols, &col_pair);
        update_mea_sets_with_cols(
          &mut mea_sets_with_cols,
          col_pair.0,
          &meas,
          &right_bp_cols_with_cols,
        );
      }
      None => {}
    }
  }
  let pseudo_col_pair = (T::zero(), sa_len + T::one());
  let mut bp_col_pairs = SparsePosMat::<T>::default();
  traceback_alifold(&mut bp_col_pairs, &pseudo_col_pair, &mea_sets_with_cols);
  bp_col_pairs
}

pub fn traceback_alifold<T>(
  bp_col_pairs: &mut SparsePosMat<T>,
  col_pair: &PosPair<T>,
  mea_sets_with_cols: &MeaSetsWithPoss<T>,
) where
  T: HashIndex,
{
  let mut mea;
  let meas = get_meas(&mea_sets_with_cols, &col_pair);
  let (i, j) = *col_pair;
  let mut k = j - T::one();
  while k > i {
    mea = meas[&k];
    let ea = meas[&(k - T::one())];
    if ea == mea {
      k = k - T::one();
      continue;
    }
    match mea_sets_with_cols.get(&k) {
      Some(meas_4_bps) => {
        for (&col_left, mea_4_bp) in meas_4_bps {
          if !(i < col_left) {
            continue;
          }
          let col_4_bp = col_left - T::one();
          let ea = meas[&col_4_bp];
          let ea = ea + mea_4_bp;
          if ea == mea {
            let col_pair = (col_left, k);
            traceback_alifold(bp_col_pairs, &col_pair, mea_sets_with_cols);
            let col_pair = (col_left - T::one(), k - T::one());
            bp_col_pairs.insert(col_pair);
            k = col_4_bp;
            break;
          }
        }
      }
      None => {}
    }
  }
}

pub fn update_mea_sets_with_cols<T>(
  mea_sets_with_cols: &mut MeaSetsWithPoss<T>,
  i: T,
  meas: &SparseProbs<T>,
  right_bp_cols_with_cols: &ColSetsWithCols<T>,
) where
  T: HashIndex,
{
  let ref right_bp_cols = right_bp_cols_with_cols[&i];
  for (&j, &weight) in right_bp_cols {
    let mea_4_bp = weight + meas[&(j - T::one())];
    match mea_sets_with_cols.get_mut(&j) {
      Some(meas_4_bps) => {
        meas_4_bps.insert(i, mea_4_bp);
      }
      None => {
        let mut meas_4_bps = SparseProbs::default();
        meas_4_bps.insert(i, mea_4_bp);
        mea_sets_with_cols.insert(j, meas_4_bps);
      }
    }
  }
}

pub fn get_meas<T>(mea_sets_with_cols: &MeaSetsWithPoss<T>, col_pair: &PosPair<T>) -> SparseProbs<T>
where
  T: HashIndex,
{
  let (i, j) = *col_pair;
  let mut meas = SparseProbs::<T>::default();
  for k in range(i, j) {
    if k == i {
      meas.insert(k, 0.);
      continue;
    }
    let mut mea = meas[&(k - T::one())];
    match mea_sets_with_cols.get(&k) {
      Some(meas_4_bps) => {
        for (&l, mea_4_bp) in meas_4_bps {
          if !(i < l) {
            continue;
          }
          let ea = meas[&(l - T::one())];
          let ea = ea + mea_4_bp;
          if ea > mea {
            mea = ea;
          }
        }
      }
      None => {}
    }
    meas.insert(k, mea);
  }
  meas
}

pub fn get_bpp_mat_alifold<T>(sa_file_path: &Path) -> SparseProbMat<T>
where
  T: HashIndex,
{
  let sa_file_prefix = sa_file_path.file_stem().unwrap().to_str().unwrap();
  let arg = format!("--id-prefix={}", sa_file_prefix);
  let args = vec![
    "-p",
    sa_file_path.to_str().unwrap(),
    &arg,
    "--noPS",
    "--noDP",
  ];
  let _ = run_command("RNAalifold", &args, "Failed to run RNAalifold");
  let mut bpp_mat_alifold = SparseProbMat::<T>::default();
  let cwd = env::current_dir().unwrap();
  let output_file_path = cwd.join(String::from(sa_file_prefix) + "_0001_ali.out");
  let output_file = BufReader::new(File::open(output_file_path.clone()).unwrap());
  for (k, line) in output_file.lines().enumerate() {
    if k == 0 {
      continue;
    }
    let line = line.unwrap();
    if !line.starts_with(" ") {
      continue;
    }
    let substrings: Vec<&str> = line.split_whitespace().collect();
    let i = T::from_usize(substrings[0].parse().unwrap()).unwrap() - T::one();
    let j = T::from_usize(substrings[1].parse().unwrap()).unwrap() - T::one();
    let mut bpp = String::from(substrings[3]);
    bpp.pop();
    let bpp = 0.01 * bpp.parse::<Prob>().unwrap();
    if bpp == 0. {
      continue;
    }
    bpp_mat_alifold.insert((i, j), bpp);
  }
  let _ = remove_file(output_file_path);
  bpp_mat_alifold
}

pub fn get_mix_bpp_mat<T>(
  sa: &SeqAlign<T>,
  bpp_mats: &SparseProbMats<T>,
  bpp_mat_alifold: &SparseProbMat<T>,
  mix_weight: Prob,
) -> SparseProbMat<T>
where
  T: HashIndex,
{
  let mut mix_bpp_mat = SparseProbMat::<T>::default();
  let sa_len = sa.cols.len();
  let num_of_rnas = sa.cols[0].len();
  for i in 0..sa_len {
    let ref pos_maps = sa.pos_map_sets[i];
    let short_i = T::from_usize(i).unwrap();
    for j in i + 1..sa_len {
      let short_j = T::from_usize(j).unwrap();
      let pos_pair = (short_i, short_j);
      let bpp_alifold = match bpp_mat_alifold.get(&pos_pair) {
        Some(&bpp_alifold) => bpp_alifold,
        None => 0.,
      };
      let ref pos_maps_2 = sa.pos_map_sets[j];
      let pos_map_pairs: Vec<(T, T)> = pos_maps
        .iter()
        .zip(pos_maps_2.iter())
        .map(|(&x, &y)| (x, y))
        .collect();
      let mut bpp_sum = 0.;
      for (pos_map_pair, bpp_mat) in pos_map_pairs.iter().zip(bpp_mats.iter()) {
        match bpp_mat.get(pos_map_pair) {
          Some(&bpp) => {
            bpp_sum += bpp;
          }
          None => {}
        }
      }
      let bpp_avg = bpp_sum / num_of_rnas as Prob;
      let mix_bpp = mix_weight * bpp_avg + (1. - mix_weight) * bpp_alifold;
      let pos_pair = (pos_pair.0 + T::one(), pos_pair.1 + T::one());
      mix_bpp_mat.insert(pos_pair, mix_bpp);
    }
  }
  mix_bpp_mat
}

pub fn revert_char(c: Base) -> u8 {
  match c {
    A => BIG_A,
    C => BIG_C,
    G => BIG_G,
    U => BIG_U,
    PSEUDO_BASE => GAP,
    _ => {
      assert!(false);
      GAP
    }
  }
}

pub fn run_command(command: &str, args: &[&str], expect: &str) -> Output {
  Command::new(command).args(args).output().expect(expect)
}
