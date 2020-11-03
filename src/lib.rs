extern crate consprob;
extern crate conshomfold;

pub use consprob::*;
pub use conshomfold::*;
pub use conshomfold::{RnaId, RnaIdPair, Prob4dMat, SparseProbMat};

// pub type Prob4dMatsWithRnaIdPairs = HashMap<RnaIdPair, Prob4dMat>;
pub type Prob4dMatsWithRnaIdPairs<T> = HashMap<RnaIdPair, Prob4dMat<T>>;
pub type ProbsWithRnaIds = Vec<Probs>;
// pub type ProbsWithPosPairs = HashMap<PosPair, Prob>;
pub type ProbsWithPosPairs<T> = HashMap<PosPair<T>, Prob>;
// pub type ProbMatsWithRnaIds = Vec<SparseProbMat>;
pub type ProbMatsWithRnaIds<T> = Vec<SparseProbMat<T>>;
pub type Col = Vec<Char>;
pub type Cols = Vec<Col>;
// pub type PosMaps = Vec<Pos>;
pub type PosMaps<T> = Vec<T>;
// pub type PosMapSets = Vec<PosMaps>;
pub type PosMapSets<T> = Vec<PosMaps<T>>;
#[derive(Debug)]
pub struct SeqAlign<T> {
  pub cols: Cols,
  pub pos_map_sets: PosMapSets<T>,
}
pub struct MeaCss<T> {
  pub bpa_pos_pairs: PosPairs<T>,
  pub ea: Mea,
}
pub type SeqId = String;
pub type SeqIds = Vec<SeqId>;

impl<T> SeqAlign<T> {
  pub fn new() -> SeqAlign<T> {
    SeqAlign {
      cols: Cols::new(),
      pos_map_sets: PosMapSets::<T>::new(),
    }
  }
}

impl<T> MeaCss<T> {
  pub fn new() -> MeaCss<T> {
    MeaCss {
      bpa_pos_pairs: PosPairs::<T>::new(),
      ea: 0.,
    }
  }
}

pub const GAP: Char = '-' as Char;

pub fn consalifold<T>(mix_bpp_mat: &ProbMat, gamma: Prob, sa: &SeqAlign<T>) -> MeaCss<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
// pub fn consalifold(mix_bpp_mat: &ProbMat, gamma: Prob, sa: &SeqAlign) -> MeaCss {
  let sa_len = sa.cols.len();
  let mut mea_mat = vec![vec![0.; sa_len]; sa_len];
  // let sa_len = sa_len as Pos;
  let sa_len = T::from_usize(sa_len).unwrap();
  let gamma_plus_1 = gamma + 1.;
  // for sub_sa_len in 1 .. sa_len + 1 {
  for sub_sa_len in range_inclusive(T::one(), sa_len) {
    // for i in 0 .. sa_len + 1 - sub_sa_len {
    for i in range_inclusive(T::zero(), sa_len - sub_sa_len) {
      // let j = i + sub_sa_len - 1;
      let j = i + sub_sa_len - T::one();
      // let (long_i, long_j) = (i as usize, j as usize);
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      if i == j {
        continue;
      }
      let mut mea = mea_mat[long_i + 1][long_j];
      let ea = mea_mat[long_i][long_j - 1];
      if ea > mea {
        mea = ea;
      }
      let ea = mea_mat[long_i + 1][long_j - 1] + gamma_plus_1 * mix_bpp_mat[long_i][long_j] - 1.;
      if ea > mea {
        mea = ea;
      }
      for k in long_i .. long_j {
        let ea = mea_mat[long_i][k] + mea_mat[k + 1][long_j];
        if ea > mea {
          mea = ea;
        }
      }
      mea_mat[long_i][long_j] = mea;
    }
  }
  let mut mea_css = MeaCss::new();
  // let mut pos_pair_stack = vec![(0, sa_len - 1)];
  let mut pos_pair_stack = vec![(T::zero(), sa_len - T::one())];
  while pos_pair_stack.len() > 0 {
    let pos_pair = pos_pair_stack.pop().expect("Failed to pop an element of a vector.");
    let (i, j) = pos_pair;
    if j <= i {continue;}
    // let (long_i, long_j) = (i as usize, j as usize);
    let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
    let mea = mea_mat[long_i][long_j];
    if mea == mea_mat[long_i + 1][long_j] {
      // pos_pair_stack.push((i + 1, j));
      pos_pair_stack.push((i + T::one(), j));
    } else if mea == mea_mat[long_i][long_j - 1] {
      // pos_pair_stack.push((i, j - 1));
      pos_pair_stack.push((i, j - T::one()));
    } else if mea == mea_mat[long_i + 1][long_j - 1] + gamma_plus_1 * mix_bpp_mat[long_i][long_j] - 1. {
      // pos_pair_stack.push((i + 1, j - 1));
      pos_pair_stack.push((i + T::one(), j - T::one()));
      mea_css.bpa_pos_pairs.push(pos_pair);
    } else {
      // for k in i .. j {
      for k in range(i, j) {
        // let long_k = k as usize;
        let long_k = k.to_usize().unwrap();
        if mea == mea_mat[long_i][long_k] + mea_mat[long_k + 1][long_j] {
          pos_pair_stack.push((i, k));
          // pos_pair_stack.push((k + 1, j));
          pos_pair_stack.push((k + T::one(), j));
          break;
        }
      }
    }
  }
  // mea_css.ea = mea_mat[0][sa_len as usize - 1];
  mea_css.ea = mea_mat[0][sa_len.to_usize().unwrap() - 1];
  mea_css
}
