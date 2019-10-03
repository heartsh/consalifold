extern crate rnafamprob;
extern crate neofold;
/* extern crate petgraph;
extern crate itertools; */

pub use rnafamprob::*;
pub use neofold::*;
/* pub use petgraph::graph::{Graph, NodeIndex};
pub use petgraph::Directed;
use itertools::Itertools; */

type MeaMat = Vec<Meas>;
// pub type Mea = Prob;
/* pub type Col = Pos;
pub type ColQuadruple = (Col, Col, Col, Col);
pub type ColQuadruples = Vec<ColQuadruple>; 
pub type ColQuadrupleSeqsWithColQuadruples = HashMap<ColQuadruple, ColQuadruples, Hasher>; 
pub type ColPair = (Col, Col);
pub type ColPairs = Vec<ColPair>;
pub type ColPairSeqsWithColPairs = HashMap<ColPair, ColPairs, Hasher>;
pub type RnaIdPosTriple = (RnaId, Pos, Pos);
pub type RnaIdPosTriples = Vec<RnaIdPosTriple>;
pub type RnaIdPosTripleSeqsWithColPairs = HashMap<ColPair, RnaIdPosTriples, Hasher>;
pub type RnaIds = Vec<RnaId>;
#[derive(Clone)]
pub struct MeaCss {
  // pub corresponding_col_quadruple_seqs_inside_col_quadruples: ColQuadrupleSeqsWithColQuadruples,
  pub seq_num: usize,
  pub mea: Mea,
  // pub rna_id_pos_pair_seqs_with_cols: RnaIdPosPairSeqs,
  // pub rna_id_pos_triple_seqs_with_col_pairs: RnaIdPosTripleSeqsWithColPairs,
  // pub pos_pair_seqs_with_col_pairs: FloatPosPairSeqsWithColPairs,
  pub pos_pair_seqs: FloatPosPairSeqs,
  pub pos_seqs_with_cols: FloatPosSeqs,
  pub col_num: usize,
  pub rna_ids: RnaIds,
  // pub pseudo_col_quadruple: ColQuadruple,
} */
pub type Prob4dMatsWithRnaIdPairs = HashMap<RnaIdPair, Prob4dMat, Hasher>;
pub type ProbsWithRnaIds = Vec<Probs>;
pub type ProbsWithPosPairs = HashMap<PosPair, Prob, Hasher>;
/* pub type MeaCssPair<'a> = (&'a MeaCss, &'a MeaCss);
pub type Meas = Vec<Mea>;
pub type MeaMat = Vec<Meas>;
pub type Mea4dMat = HashMap<ColQuadruple, Mea, Hasher>;
pub type ClusterIndex = RnaId;
pub type ClusterIndexes = Vec<ClusterIndex>;
pub type ClusterScore = Mea;
pub type GuideTree = Graph<ClusterScore, ClusterScore, Directed, usize>;
type ClusterIndexPair = (ClusterIndex, ClusterIndex);
type ClusterScoreMat = HashMap<ClusterIndexPair, ClusterScore, Hasher>;
struct ClusterIndexScorePair {
  pub index: ClusterIndex,
  pub score: ClusterScore,
}
type ClusterIndexScorePairs = Vec<ClusterIndexScorePair>;
type ClusterIndexScorePairSeqsWithClusterIndexes = HashMap<ClusterIndex, ClusterIndexScorePairs, Hasher>;
pub type SparseMeaMat = HashMap<RnaIdPair, Mea>;
type SeqNumsWithClusterIndexes = HashMap<ClusterIndex, usize, Hasher>;
type RnaIdPosPair = (RnaId, Pos);
type RnaIdPosPairs = Vec<RnaIdPosPair>;
type RnaIdPosPairSeqs = Vec<RnaIdPosPairs>;
pub type BoolsWithPosPairs = HashMap<PosPair, bool, Hasher>;
pub type ProbSeqsWithRnaIds = Vec<Probs>;
pub type FloatPos = f32;
pub type FloatPosPair = (FloatPos, FloatPos);
pub type FloatPosPairs = Vec<FloatPosPair>;
pub type FloatPosPairSeqsWithColPairs = HashMap<ColPair, FloatPosPairs, Hasher>;
pub type FloatPoss = Vec<FloatPos>;
pub type FloatPosSeqs = Vec<FloatPoss>;
// pub type FloatPosSeqsWithCols = HashMap<Col, FloatPoss, Hasher>;
pub type FloatPosPairSeqs = Vec<FloatPosPairs>;
pub type ProbsWithCols = HashMap<Col, Prob, Hasher>; */
pub type SeqId = String;
pub type SeqIds = Vec<SeqId>;
// type FastaRecord = (FastaId, Seq, usize);
// type FastaRecords = Vec<FastaRecord>;
pub type Col = Vec<Char>;
pub type Cols = Vec<Col>;
pub type PosMaps = Vec<Pos>;
pub type PosMapSets = Vec<PosMaps>;
#[derive(Debug)]
pub struct SeqAlign {
  // pub seq_ids: SeqIds,
  pub cols: Cols,
  pub pos_map_sets: PosMapSets,
}
pub struct MeaCss {
  pub bpa_pos_pairs: PosPairs,
  pub ea: Mea,
}

impl SeqAlign {
  pub fn new() -> SeqAlign {
    SeqAlign {
      // seq_ids: SeqIds::new(),
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

/* impl MeaCss {
  pub fn new() -> MeaCss {
    MeaCss {
      // corresponding_col_quadruple_seqs_inside_col_quadruples: ColQuadrupleSeqsWithColQuadruples::default(),
      seq_num: 0,
      mea: 0.,
      /* rna_id_pos_pair_seqs_with_cols: RnaIdPosPairSeqs::new(),
      rna_id_pos_triple_seqs_with_col_pairs: RnaIdPosTripleSeqsWithColPairs::default(), */
      pos_pair_seqs: FloatPosPairSeqs::new(),
      pos_seqs_with_cols: FloatPosSeqs::new(),
      col_num: 0,
      rna_ids: RnaIds::new(),
      // pseudo_col_quadruple: (0, 0, 0, 0),
    }
  }
} */

/* impl ClusterIndexScorePair {
  pub fn new(cluster_index: ClusterIndex, cluster_score: ClusterScore) -> ClusterIndexScorePair {
    ClusterIndexScorePair {
      index: cluster_index,
      score: cluster_score,
    }
  }
} */

pub const GAP: Char = '-' as Char;

#[inline]
pub fn neoalifold(bpap_mats_with_rna_id_pairs: &Prob4dMatsWithRnaIdPairs, mean_upp_mat: &Probs, gamma: Prob, sa: &SeqAlign) -> MeaCss {
  let sa_len = sa.cols.len();
  // let mut mea_mat_4_bpa_pos_pairs = MeaMat::default();
  let mut mea_mat = vec![vec![0.; sa_len]; sa_len];
  // let mut pos_seqs_with_poss_4_forward_bpas = PosSeqsWithPoss::default();
  let num_of_rnas = sa.cols[0].len();
  let combination_num = (num_of_rnas * (num_of_rnas - 1)) as Prob / 2.;
  let mut mean_bpaps_with_pos_pairs = ProbsWithPosPairs::default();
  for sub_sa_len in 2 .. sa_len + 1 {
    for i in 0 .. sa_len + 1 - sub_sa_len {
      let j = i + sub_sa_len - 1;
      let pos_pair = (i, j);
      let mut mea = mea_mat[i + 1][j] + mean_upp_mat[i];
      // println!("MEA 1: {}", mea);
      let ea = mea_mat[i][j - 1] + mean_upp_mat[j];
      if ea > mea {
        mea = ea;
      }
      // println!("MEA 2: {}", mea);
      let mut mean_bpap = 0.;
      for rna_id_1 in 0 .. num_of_rnas {
        for rna_id_2 in rna_id_1 + 1 .. num_of_rnas {
          let rna_id_pair = (rna_id_1, rna_id_2);
          let ref bpap_mat = bpap_mats_with_rna_id_pairs[&rna_id_pair];
          let base_quadruple = (sa.cols[i][rna_id_1], sa.cols[j][rna_id_1], sa.cols[i][rna_id_2], sa.cols[j][rna_id_2]);
          if base_quadruple.0 == GAP || base_quadruple.1 == GAP || base_quadruple.2 == GAP || base_quadruple.3 == GAP {
            continue;
          }
          let pos_quadruple = (sa.pos_map_sets[i][rna_id_1], sa.pos_map_sets[j][rna_id_1], sa.pos_map_sets[i][rna_id_2], sa.pos_map_sets[j][rna_id_2]);
          match bpap_mat.get(&pos_quadruple) {
            Some(&bpap) => {
              mean_bpap += bpap;
            },
            None => {},
          }
        }
      }
      mean_bpap /= combination_num;
      mean_bpaps_with_pos_pairs.insert(pos_pair, mean_bpap);
      let ea = mea_mat[i + 1][j - 1] + gamma * mean_bpap;
      if ea > mea {
        mea = ea;
      }
      // println!("MEA 3: {}", mea);
      for k in i .. j {
        let ea = mea_mat[i][k] + mea_mat[k + 1][j];
        if ea > mea {
          mea = ea;
        }
        // println!("MEA 4: {}", mea);
      }
      mea_mat[i][j] = mea;
    }
  }
  let mut mea_css = MeaCss::new();
  let pseudo_pos_pair = (0, sa_len - 1);
  let mut pos_pair_stack = vec![pseudo_pos_pair];
  while pos_pair_stack.len() > 0 {
    let pos_pair = pos_pair_stack.pop().expect("Failed to pop an element of a vector.");
    let (i, j) = pos_pair;
    let mea = mea_mat[i][j];
    if mea == 0. {continue;}
    let contains_mean_bpap = mean_bpaps_with_pos_pairs.contains_key(&pos_pair);
    let mean_bpap = if contains_mean_bpap {mean_bpaps_with_pos_pairs[&pos_pair]} else {0.};
    if mea == mea_mat[i + 1][j] + mean_upp_mat[i] {
      pos_pair_stack.push((i + 1, j));
    } else if mea == mea_mat[i][j - 1] + mean_upp_mat[j] {
      pos_pair_stack.push((i, j - 1));
    } else if contains_mean_bpap && mea == mea_mat[i + 1][j - 1] + gamma * mean_bpap {
      pos_pair_stack.push((i + 1, j - 1));
      mea_css.bpa_pos_pairs.push(pos_pair);
    } else {
      for k in i .. j {
        if mea == mea_mat[i][k] + mea_mat[k + 1][j] {
          pos_pair_stack.push((i, k));
          pos_pair_stack.push((k + 1, j));
          break;
        }
      }
    }
  }
  mea_css.ea = mea_mat[0][sa_len - 1];
  mea_css
}

/* #[inline]
pub fn get_guide_tree(mea_mat: &SparseMeaMat, seq_num: usize) -> GuideTree {
  let mut guide_tree = GuideTree::default();
  for _ in 0 .. seq_num {
    guide_tree.add_node(NEG_INFINITY);
  }
  let mut cluster_score_mat = ClusterScoreMat::default();
  let mut cluster_index_score_pair_seqs_with_cluster_indexes = ClusterIndexScorePairSeqsWithClusterIndexes::default();
  for i in 0 .. seq_num {
    let mut cluster_index_score_pairs = ClusterIndexScorePairs::new();
    for j in i + 1 .. seq_num {
      let mea = mea_mat[&(i, j)];
      cluster_score_mat.insert((i, j), mea);
      cluster_index_score_pairs.push(ClusterIndexScorePair::new(j, mea));
    }
    cluster_index_score_pair_seqs_with_cluster_indexes.insert(i, cluster_index_score_pairs);
  }
  for cluster_index_score_pairs in cluster_index_score_pair_seqs_with_cluster_indexes.values_mut() {
    cluster_index_score_pairs.sort_unstable_by(|cluster_index_score_pair_1, cluster_index_score_pair_2| {cluster_index_score_pair_2.score.partial_cmp(&cluster_index_score_pair_1.score).expect("Failed to compare 2 floating-point numbers with each other.")});
  }
  let mut seq_nums_with_cluster_indexes = (0 .. seq_num).map(|cluster_index| {(cluster_index, 1)}).collect::<SeqNumsWithClusterIndexes>();
  let mut new_cluster_index = seq_num;
  while seq_nums_with_cluster_indexes.len() > 1 {
    let mut cluster_index_pair_4_merge = (0, 0);
    let mut max_cluster_score = NEG_INFINITY;
    for (&i, cluster_index_score_pairs) in cluster_index_score_pair_seqs_with_cluster_indexes.iter() {
      if cluster_index_score_pairs.len() == 0 {continue;}
      let ref cluster_index_score_pair = cluster_index_score_pairs[0];
      if cluster_index_score_pair.score > max_cluster_score {
        cluster_index_pair_4_merge.0 = i;
        cluster_index_pair_4_merge.1 = cluster_index_score_pair.index;
        max_cluster_score = cluster_index_score_pair.score;
      }
    }
    let pair_of_nums_of_seqs_in_clusters = (seq_nums_with_cluster_indexes[&cluster_index_pair_4_merge.0], seq_nums_with_cluster_indexes[&cluster_index_pair_4_merge.1]);
    let sum_of_seq_num_pair = pair_of_nums_of_seqs_in_clusters.0 + pair_of_nums_of_seqs_in_clusters.1;
    seq_nums_with_cluster_indexes.remove(&cluster_index_pair_4_merge.0);
    seq_nums_with_cluster_indexes.remove(&cluster_index_pair_4_merge.1);
    for &cluster_index in seq_nums_with_cluster_indexes.keys() {
      let ordered_cluster_index_pair_1 = (min(cluster_index_pair_4_merge.0, cluster_index), max(cluster_index_pair_4_merge.0, cluster_index));
      let ordered_cluster_index_pair_2 = (min(cluster_index_pair_4_merge.1, cluster_index), max(cluster_index_pair_4_merge.1, cluster_index));
      let new_cluster_score = (pair_of_nums_of_seqs_in_clusters.0 as ClusterScore * cluster_score_mat[&ordered_cluster_index_pair_1]
      + pair_of_nums_of_seqs_in_clusters.1 as ClusterScore * cluster_score_mat[&ordered_cluster_index_pair_2])
      / sum_of_seq_num_pair as ClusterScore;
      cluster_score_mat.remove(&ordered_cluster_index_pair_1);
      cluster_score_mat.remove(&ordered_cluster_index_pair_2);
      cluster_score_mat.insert((cluster_index, new_cluster_index), new_cluster_score);
      let cluster_index_score_pairs = cluster_index_score_pair_seqs_with_cluster_indexes.get_mut(&cluster_index).expect("Failed to get an element from a hash map.");
      let mut indexs_4_remove = Vec::new();
      for (i, cluster_index_score_pair) in cluster_index_score_pairs.iter().enumerate() {
        if cluster_index_score_pair.index == cluster_index_pair_4_merge.0 || cluster_index_score_pair.index == cluster_index_pair_4_merge.1 {
          indexs_4_remove.push(i);
        }
      }
      if indexs_4_remove.len() > 0 {
        cluster_index_score_pairs.remove(indexs_4_remove[0]);
      }
      if indexs_4_remove.len() > 1 {
        cluster_index_score_pairs.remove(indexs_4_remove[1] - 1);
      }
      cluster_index_score_pairs.push(ClusterIndexScorePair::new(new_cluster_index, new_cluster_score));
      cluster_index_score_pairs.sort_unstable_by(|cluster_index_score_pair_1, cluster_index_score_pair_2| {cluster_index_score_pair_2.score.partial_cmp(&cluster_index_score_pair_1.score).expect("Failed to compare 2 floating-point numbers.")});
    }
    cluster_score_mat.remove(&cluster_index_pair_4_merge);
    cluster_index_score_pair_seqs_with_cluster_indexes.remove(&cluster_index_pair_4_merge.0);
    cluster_index_score_pair_seqs_with_cluster_indexes.remove(&cluster_index_pair_4_merge.1);
    cluster_index_score_pair_seqs_with_cluster_indexes.insert(new_cluster_index, ClusterIndexScorePairs::new());
    seq_nums_with_cluster_indexes.insert(new_cluster_index, sum_of_seq_num_pair);
    let new_node = guide_tree.add_node(max_cluster_score);
    guide_tree.add_edge(new_node, NodeIndex::new(cluster_index_pair_4_merge.0), max_cluster_score);
    guide_tree.add_edge(new_node, NodeIndex::new(cluster_index_pair_4_merge.1), max_cluster_score);
    new_cluster_index += 1;
  }
  guide_tree
}

#[inline]
pub fn is_gap_pos(pos: FloatPos) -> bool {
  if pos.fract() != 0. {true} else {false}
} */
