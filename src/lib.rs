extern crate stem;
extern crate petgraph;
extern crate itertools;

pub use stem::*;
pub use petgraph::graph::{Graph, NodeIndex};
pub use petgraph::Directed;
use itertools::Itertools;

pub type Mea = Prob;
pub type Col = Pos;
pub type ColQuadruple = (Col, Col, Col, Col);
pub type ColQuadruples = Vec<ColQuadruple>; 
pub type ColQuadrupleSeqsWithColQuadruples = HashMap<ColQuadruple, ColQuadruples, Hasher>; 
pub type ColPair = (Col, Col);
pub type ColPairs = Vec<ColPair>;
pub type ColPairSeqsWithColPairs = HashMap<ColPair, ColPairs, Hasher>; 
pub type ProbsWithColPairs = HashMap<ColPair, Prob, Hasher>;
pub type PosPairs = Vec<PosPair>;
pub type PosPairSeqsWithColPairs = HashMap<ColPair, PosPairs, Hasher>;
pub type RnaIds = Vec<RnaId>;
#[derive(Clone)]
pub struct MeaCss {
  pub corresponding_col_quadruple_seqs_inside_col_quadruples: ColQuadrupleSeqsWithColQuadruples,
  pub mean_bpaps_with_col_pairs: ProbsWithColPairs,
  pub seq_num: usize,
  pub mea: Mea,
  pub pos_pair_seqs_with_col_pairs: PosPairSeqsWithColPairs,
  pub col_num: usize,
  pub rna_ids: RnaIds,
  pub pseudo_col_quadruple: ColQuadruple,
}
pub type Prob4dMatsWithRnaIdPairs = HashMap<RnaIdPair, Prob4dMat, Hasher>;
pub type MeaCssPair<'a> = (&'a MeaCss, &'a MeaCss);
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
type Poss = Vec<Pos>;
type PosSeqs = Vec<Poss>;
type BoolsWithPosPairs = HashMap<PosPair, bool, Hasher>;

impl MeaCss {
  pub fn new() -> MeaCss {
    MeaCss {
      corresponding_col_quadruple_seqs_inside_col_quadruples: ColQuadrupleSeqsWithColQuadruples::default(),
      mean_bpaps_with_col_pairs: ProbsWithColPairs::default(),
      seq_num: 0,
      mea: 0.,
      pos_pair_seqs_with_col_pairs: PosPairSeqsWithColPairs::default(),
      col_num: 0,
      rna_ids: RnaIds::new(),
      pseudo_col_quadruple: (0, 0, 0, 0),
    }
  }
}

impl ClusterIndexScorePair {
  pub fn new(cluster_index: ClusterIndex, cluster_score: ClusterScore) -> ClusterIndexScorePair {
    ClusterIndexScorePair {
      index: cluster_index,
      score: cluster_score,
    }
  }
}

#[inline]
pub fn get_mea_consensus_ss(mea_css_pair: &MeaCssPair, gamma_plus_1: Prob, bpap_mats_with_rna_id_pairs: &Prob4dMatsWithRnaIdPairs) -> MeaCss {
  let mut mea_mat_4_corresponding_col_quadruples = Mea4dMat::default();
  let mut col_pair_seqs_with_col_pairs_4_forward_bpas = ColPairSeqsWithColPairs::default();
  let inverse_gamma_plus_1 = 1. / gamma_plus_1;
  let seq_num_pair = (mea_css_pair.0.seq_num, mea_css_pair.1.seq_num);
  let combi_num_pair = (
    (0 .. seq_num_pair.0).combinations(2).fold(0, |acc, i| {&acc + 1}) as Prob,
    (0 .. seq_num_pair.1).combinations(2).fold(0, |acc, i| {&acc + 1}) as Prob,
  );
  let sum_of_seq_num_pair = seq_num_pair.0 + seq_num_pair.1;
  for sub_seq_len_1 in 2 .. mea_css_pair.0.col_num + 1 {
    for i in 0 .. mea_css_pair.0.col_num - sub_seq_len_1 + 1 {
      let j = i + sub_seq_len_1 - 1;
      if !(mea_css_pair.0.seq_num == 1 || mea_css_pair.0.mean_bpaps_with_col_pairs.contains_key(&(i, j))) {continue;}
      for sub_seq_len_2 in 2 .. mea_css_pair.1.col_num + 1 {
        for k in 0 .. mea_css_pair.1.col_num - sub_seq_len_2 + 1 {
          let l = k + sub_seq_len_2 - 1;
          let col_quadruple = (i, j, k, l);
          if !(mea_css_pair.1.seq_num == 1 || mea_css_pair.1.mean_bpaps_with_col_pairs.contains_key(&(k, l))) {continue;}
          let mut mean_bpap = if mea_css_pair.0.seq_num == 1 {0.} else {combi_num_pair.0 * mea_css_pair.0.mean_bpaps_with_col_pairs[&(i, j)]} + if mea_css_pair.1.seq_num == 1 {0.} else {combi_num_pair.1 * mea_css_pair.1.mean_bpaps_with_col_pairs[&(k, l)]};
          let ref pos_pairs_1 = mea_css_pair.0.pos_pair_seqs_with_col_pairs[&(i, j)];
          let ref pos_pairs_2 = mea_css_pair.1.pos_pair_seqs_with_col_pairs[&(k, l)];
          for (m, pos_pair_1) in pos_pairs_1.iter().enumerate() {
            let rna_id_1 = mea_css_pair.0.rna_ids[m];
            for (n, pos_pair_2) in pos_pairs_2.iter().enumerate() {
              let rna_id_2 = mea_css_pair.1.rna_ids[n];
              let rna_id_pair = if rna_id_1 < rna_id_2 {(rna_id_1, rna_id_2)} else {(rna_id_2, rna_id_1)};
              let pos_quadruple = if rna_id_1 < rna_id_2 {(pos_pair_1.0, pos_pair_1.1, pos_pair_2.0, pos_pair_2.1)} else {(pos_pair_2.0, pos_pair_2.1, pos_pair_1.0, pos_pair_1.1)};
              let ref bpap_mat = bpap_mats_with_rna_id_pairs[&rna_id_pair];
              if !bpap_mat.contains_key(&pos_quadruple) {continue;}
              mean_bpap += bpap_mat[&pos_quadruple];
            }
          }
          mean_bpap /= (0 .. sum_of_seq_num_pair).combinations(2).fold(0, |acc, i| {&acc + 1}) as Prob;
          if mean_bpap <= inverse_gamma_plus_1 {
            continue;
          }
          let mea_mat_4_corresponding_col_quadruple = get_mea_mat_4_corresponding_col_quadruple(&col_quadruple, &mea_mat_4_corresponding_col_quadruples, &col_pair_seqs_with_col_pairs_4_forward_bpas);
          mea_mat_4_corresponding_col_quadruples.insert(col_quadruple, mea_mat_4_corresponding_col_quadruple[j - i - 1][l - k - 1] + gamma_plus_1 * mean_bpap - 1.);
          if col_pair_seqs_with_col_pairs_4_forward_bpas.contains_key(&(j, l)) {
            col_pair_seqs_with_col_pairs_4_forward_bpas.get_mut(&(j, l)).expect("Failed to get an element of a hash map.").push((i, k));
          } else {
            col_pair_seqs_with_col_pairs_4_forward_bpas.insert((j, l), vec![(i, k)]);
          }
        }
      }
    }
  }
  let mut mea_css = MeaCss::new();
  let pseudo_col_quadruple = (0, mea_css_pair.0.col_num - 1, 0, mea_css_pair.1.col_num - 1);
  let mut col_quadruple_stack = vec![pseudo_col_quadruple];
  let mut pos_seqs = PosSeqs::new();
  let mut pos_pairs_exist_in_fisrt_seq = BoolsWithPosPairs::default();
  while col_quadruple_stack.len() > 0 {
    let col_quadruple_1 = col_quadruple_stack.pop().expect("Failed to pop an element of a vector.");
    let (i, j, k, l) = col_quadruple_1;
    pos_pairs_exist_in_fisrt_seq.insert(mea_css_pair.0.pos_pair_seqs_with_col_pairs[&(i, j)][0], true);
    let mut left_poss = Poss::new();
    let mut right_poss = Poss::new();
    for &(m, n) in &mea_css_pair.0.pos_pair_seqs_with_col_pairs[&(i, j)] {
      left_poss.push(m);
      right_poss.push(n);
    }
    for &(m, n) in &mea_css_pair.1.pos_pair_seqs_with_col_pairs[&(k, l)] {
      left_poss.push(m);
      right_poss.push(n);
    }
    pos_seqs.push(left_poss);
    pos_seqs.push(right_poss);
    mea_css.col_num += 2;
    let mea_mat_4_corresponding_col_quadruple = get_mea_mat_4_corresponding_col_quadruple(&col_quadruple_1, &mea_mat_4_corresponding_col_quadruples, &col_pair_seqs_with_col_pairs_4_forward_bpas);
    let mea = mea_mat_4_corresponding_col_quadruple[j - i - 1][l - k - 1];
    if mea == 0. {continue;}
    let mut n = j - 1;
    let mut p = l - 1;
    while mea_mat_4_corresponding_col_quadruple[n - i][p - k] > 0. {
      let mea = mea_mat_4_corresponding_col_quadruple[n - i][p - k];
      if mea == mea_mat_4_corresponding_col_quadruple[n - i - 1][p - k] {
        n = n - 1;
      } else if mea == mea_mat_4_corresponding_col_quadruple[n - i][p - k - 1] {
        p = p - 1;
      } else {
        match col_pair_seqs_with_col_pairs_4_forward_bpas.get(&(n, p)) {
          Some(col_pairs) => {
            for &(m, o) in col_pairs {
              if m <= i || o <= k {continue;}
              let col_quadruple_2 = (m, n, o, p);
              if mea == mea_mat_4_corresponding_col_quadruple[m - i - 1][o - k - 1] + mea_mat_4_corresponding_col_quadruples[&col_quadruple_2] {
                if mea_css.corresponding_col_quadruple_seqs_inside_col_quadruples.contains_key(&col_quadruple_1) {
                  mea_css.corresponding_col_quadruple_seqs_inside_col_quadruples.get_mut(&col_quadruple_1).expect("Failed to get an element of a hash map.").push(col_quadruple_2);
                } else {
                  mea_css.corresponding_col_quadruple_seqs_inside_col_quadruples.insert(col_quadruple_1, vec![col_quadruple_2]);
                }
                col_quadruple_stack.push(col_quadruple_2);
                n = m - 1;
                p = o - 1;
                break;
              }
            }
          },
          None => {},
        }
      }
    }
  }
  mea_css.mea = mea_mat_4_corresponding_col_quadruples[&pseudo_col_quadruple];
  mea_css.seq_num = sum_of_seq_num_pair;
  mea_css.rna_ids = mea_css_pair.0.rna_ids.clone();
  mea_css.rna_ids.extend(&mea_css_pair.1.rna_ids);
  mea_css.pseudo_col_quadruple = pseudo_col_quadruple;
  pos_seqs.sort_unstable_by(|poss_1, poss_2| {poss_1[0].cmp(&poss_2[0])});
  for i in 0 .. mea_css.col_num {
    let ref left_poss = pos_seqs[i];
    let left_col_first_pos = left_poss[0];
    for j in i + 1 .. mea_css.col_num {
      let ref right_poss = pos_seqs[j];
      let right_col_first_pos = right_poss[0];
      let first_seq_pos_pair = (left_col_first_pos, right_col_first_pos);
      if pos_pairs_exist_in_fisrt_seq.contains_key(&first_seq_pos_pair) {
        let col_pair = (i, j);
        let pos_pairs = left_poss.iter().zip(right_poss).map(|(&i, &j)| {(i, j)}).collect::<PosPairs>();
        mea_css.pos_pair_seqs_with_col_pairs.insert(col_pair, pos_pairs);
      }
    }
  }
  for (col_pair, pos_pairs) in &mea_css.pos_pair_seqs_with_col_pairs {
    let mut mean_bpap = 0.;
    for (i, pos_pair_1) in pos_pairs.iter().enumerate() {
      let rna_id_1 = mea_css.rna_ids[i];
      for (j, pos_pair_2) in (&pos_pairs[i + 1 ..]).iter().enumerate() {
        let j = j + i + 1;
        let rna_id_2 = mea_css.rna_ids[j];
        let rna_id_pair = if rna_id_1 < rna_id_2 {(rna_id_1, rna_id_2)} else {(rna_id_2, rna_id_1)};
        let pos_quadruple = if rna_id_1 < rna_id_2 {
          (pos_pair_1.0, pos_pair_1.1, pos_pair_2.0, pos_pair_2.1)
        } else {
          (pos_pair_2.0, pos_pair_2.1, pos_pair_1.0, pos_pair_1.1)
        };
        let ref bpap_mat = bpap_mats_with_rna_id_pairs[&rna_id_pair];
        if !bpap_mat.contains_key(&pos_quadruple) {continue;}
        mean_bpap += bpap_mat[&pos_quadruple];
      }
    }
    mean_bpap /= (0 .. sum_of_seq_num_pair).combinations(2).fold(0, |acc, i| {&acc + 1}) as Prob;
    mea_css.mean_bpaps_with_col_pairs.insert(*col_pair, mean_bpap);
  }
  mea_css
}

#[inline]
fn get_mea_mat_4_corresponding_col_quadruple(col_quadruple: &ColQuadruple, mea_mat_4_corresponding_col_quadruples: &Mea4dMat, col_pair_seqs_with_col_pairs_4_forward_bpas: &ColPairSeqsWithColPairs) -> MeaMat {
  let &(i, j, k, l) = col_quadruple;
  let sub_seq_len_1 = j - i + 1;
  let sub_seq_len_2 = l - k + 1;
  let mut mea_mat_4_corresponding_col_quadruple = vec![vec![0.; sub_seq_len_2 - 1]; sub_seq_len_1 - 1];
  for n in i + 2 .. j {
    for p in k + 2 .. l {
      let mut mea = 0.;
      match col_pair_seqs_with_col_pairs_4_forward_bpas.get(&(n, p)) {
        Some(col_pairs) => {
          for &(m, o) in col_pairs {
            if m <= i || o <= k {continue;}
            let ea = mea_mat_4_corresponding_col_quadruple[m - i - 1][o - k - 1] + mea_mat_4_corresponding_col_quadruples[&(m, n, o, p)];
            if ea > mea {
              mea = ea;
            }
          }
        },
        None => {},
      }
      let ea = mea_mat_4_corresponding_col_quadruple[n - i - 1][p - k];
      if ea > mea {
        mea = ea;
      }
      let ea = mea_mat_4_corresponding_col_quadruple[n - i][p - k - 1];
      if ea > mea {
        mea = ea;
      }
      mea_mat_4_corresponding_col_quadruple[n - i][p - k] = mea;
    }
  }
  mea_mat_4_corresponding_col_quadruple
}

#[inline]
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
    cluster_index_score_pairs.sort_unstable_by(|cluster_index_score_pair_1, cluster_index_score_pair_2| {cluster_index_score_pair_2.score.partial_cmp(&cluster_index_score_pair_1.score).expect("Failed to compare 2 floating-point numbers.")});
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
