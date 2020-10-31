extern crate consprob;
extern crate conshomfold;

pub use consprob::*;
pub use conshomfold::*;
pub use conshomfold::{RnaId, RnaIdPair, Prob4dMat, SparseProbMat};
pub type Prob4dMatsWithRnaIdPairs = HashMap<RnaIdPair, Prob4dMat>;
pub type ProbsWithRnaIds = Vec<Probs>;
pub type ProbsWithPosPairs = HashMap<PosPair, Prob>;
pub type ProbMatsWithRnaIds = Vec<SparseProbMat>;
pub type Col = Vec<Char>;
pub type Cols = Vec<Col>;
pub type PosMaps = Vec<Pos>;
pub type PosMapSets = Vec<PosMaps>;
#[derive(Debug)]
pub struct SeqAlign {
  pub cols: Cols,
  pub pos_map_sets: PosMapSets,
}
pub struct MeaCss {
  pub bpa_pos_pairs: PosPairs,
  pub ea: Mea,
}
pub type SeqId = String;
pub type SeqIds = Vec<SeqId>;

impl SeqAlign {
  pub fn new() -> SeqAlign {
    SeqAlign {
      cols: Cols::new(),
      pos_map_sets: PosMapSets::new(),
    }
  }
}

impl MeaCss {
  pub fn new() -> MeaCss {
    MeaCss {
      bpa_pos_pairs: PosPairs::new(),
      ea: 0.,
    }
  }
}

pub const GAP: Char = '-' as Char;

pub fn consalifold(mix_bpp_mat: &ProbMat, gamma: Prob, sa: &SeqAlign) -> MeaCss {
  let sa_len = sa.cols.len();
  let mut mea_mat = vec![vec![0.; sa_len]; sa_len];
  let sa_len = sa_len as Pos;
  let gamma_plus_1 = gamma + 1.;
  for sub_sa_len in 1 .. sa_len + 1 {
    for i in 0 .. sa_len + 1 - sub_sa_len {
      let j = i + sub_sa_len - 1;
      let (long_i, long_j) = (i as usize, j as usize);
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
  let mut pos_pair_stack = vec![(0, sa_len - 1)];
  while pos_pair_stack.len() > 0 {
    let pos_pair = pos_pair_stack.pop().expect("Failed to pop an element of a vector.");
    let (i, j) = pos_pair;
    if j <= i {continue;}
    let (long_i, long_j) = (i as usize, j as usize);
    let mea = mea_mat[long_i][long_j];
    if mea == mea_mat[long_i + 1][long_j] {
      pos_pair_stack.push((i + 1, j));
    } else if mea == mea_mat[long_i][long_j - 1] {
      pos_pair_stack.push((i, j - 1));
    } else if mea == mea_mat[long_i + 1][long_j - 1] + gamma_plus_1 * mix_bpp_mat[long_i][long_j] - 1. {
      pos_pair_stack.push((i + 1, j - 1));
      mea_css.bpa_pos_pairs.push(pos_pair);
    } else {
      for k in i .. j {
        let long_k = k as usize;
        if mea == mea_mat[long_i][long_k] + mea_mat[long_k + 1][long_j] {
          pos_pair_stack.push((i, k));
          pos_pair_stack.push((k + 1, j));
          break;
        }
      }
    }
  }
  mea_css.ea = mea_mat[0][sa_len as usize - 1];
  mea_css
}
