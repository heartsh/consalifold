extern crate phyloprob;
extern crate phylofold;

pub use phyloprob::*;
pub use phylofold::*;
pub use phylofold::{RnaId, RnaIdPair, Prob4dMat, SparseProbMat};
pub type Prob4dMatsWithRnaIdPairs = FxHashMap<RnaIdPair, Prob4dMat>;
pub type ProbsWithRnaIds = Vec<Probs>;
pub type ProbsWithPosPairs = FxHashMap<PosPair, Prob>;
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

#[inline]
pub fn neoalifold(mean_bpp_mat: &ProbMat, mean_upp_mat: &Probs, centroidalifold_bpp_mat: &ProbMat, gamma: Prob, prob_weight: Prob, sa: &SeqAlign) -> MeaCss {
  let sa_len = sa.cols.len();
  let mut centroidalifold_upp_mat = vec![0.; sa_len];
  let mut mea_mat = vec![vec![0.; sa_len]; sa_len];
  let sa_len = sa_len as Pos;
  for i in 0 .. sa_len {
    let long_i = i as usize;
    let mut centroidalifold_upp = 1.;
    for j in 0 .. sa_len {
      if i == j {continue;}
      let long_j = j as usize;
      let pos_pair = if long_i < long_j {(long_i, long_j)} else {(long_j, long_i)};
      let centroidalifold_bpp = centroidalifold_bpp_mat[pos_pair.0][pos_pair.1];
      centroidalifold_upp -= centroidalifold_bpp;
    }
    centroidalifold_upp_mat[long_i] = centroidalifold_upp;
  }
  for sub_sa_len in 1 .. sa_len + 1 {
    for i in 0 .. sa_len + 1 - sub_sa_len {
      let j = i + sub_sa_len - 1;
      let (long_i, long_j) = (i as usize, j as usize);
      if i == j {
        mea_mat[long_i][long_j] = prob_weight * mean_upp_mat[long_i] + (1. - prob_weight) * centroidalifold_upp_mat[long_i];
        continue;
      }
      let mut mea = mea_mat[long_i + 1][long_j] + prob_weight * mean_upp_mat[long_i] + (1. - prob_weight) * centroidalifold_upp_mat[long_i];
      let ea = mea_mat[long_i][long_j - 1] + prob_weight * mean_upp_mat[long_j] + (1. - prob_weight) * centroidalifold_upp_mat[long_j];
      if ea > mea {
        mea = ea;
      }
      let ea = mea_mat[long_i + 1][long_j - 1] + gamma * (prob_weight * mean_bpp_mat[long_i][long_j] + (1. - prob_weight) * centroidalifold_bpp_mat[long_i][long_j]);
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
    if mea == mea_mat[long_i + 1][long_j] + prob_weight * mean_upp_mat[long_i] + (1. - prob_weight) * centroidalifold_upp_mat[long_i] {
      pos_pair_stack.push((i + 1, j));
    } else if mea == mea_mat[long_i][long_j - 1] + prob_weight * mean_upp_mat[long_j] + (1. - prob_weight) * centroidalifold_upp_mat[long_j] {
      pos_pair_stack.push((i, j - 1));
    } else if mea == mea_mat[long_i + 1][long_j - 1] + gamma * (prob_weight * mean_bpp_mat[long_i][long_j] + (1. - prob_weight) * centroidalifold_bpp_mat[long_i][long_j]) {
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
