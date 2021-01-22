extern crate consprob;

pub use consprob::*;

pub type Prob4dMatsWithRnaIdPairs<T> = HashMap<RnaIdPair, Prob4dMat<T>>;
pub type ProbsWithRnaIds = Vec<Probs>;
pub type ProbsWithPosPairs<T> = HashMap<PosPair<T>, Prob>;
pub type ProbMatsWithRnaIds<T> = Vec<SparseProbMat<T>>;
pub type Col = Vec<Base>;
pub type Cols = Vec<Col>;
pub type PosMaps<T> = Vec<T>;
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
pub struct CssPartFuncMats<T: Hash> {
  pub part_func_mat: PartFuncMat,
  pub part_func_mat_4_rightmost_base_pairings_on_el: PartFuncMat,
  pub part_func_mat_4_rightmost_base_pairings_on_mls: PartFuncMat,
  pub part_func_mat_4_base_pairings: SparsePartFuncMat<T>,
  pub part_func_mat_4_base_pairings_accessible_on_el: SparsePartFuncMat<T>,
  pub part_func_mat_4_base_pairings_accessible_on_mls: SparsePartFuncMat<T>,
  pub part_func_mat_4_base_unpairings_on_sa: PartFuncMat,
  pub part_func_mat_4_base_unpairings_on_el: PartFuncMat,
  pub part_func_mat_4_base_unpairings_on_mls: PartFuncMat,
  pub part_func_mat_4_at_least_1_base_pairings_on_mls: PartFuncMat,
  pub score_mat_4_ml_closing_basepairings: SparsePartFuncMat<T>,
  pub score_mat_4_loop_aligns: PartFuncs,
  pub twoloop_score_4d_mat: PartFunc4dMat<T>,
}

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

impl<T: Hash + Clone> CssPartFuncMats<T> {
  fn new(sa_len: usize) -> CssPartFuncMats<T> {
    let zero_mat = vec![vec![0.; sa_len]; sa_len];
    let neg_inf_mat = vec![vec![NEG_INFINITY; sa_len]; sa_len];
    let part_func_mat = SparsePartFuncMat::<T>::default();
    CssPartFuncMats {
      part_func_mat: vec![vec![0.; sa_len]; sa_len],
      part_func_mat_4_rightmost_base_pairings_on_el: neg_inf_mat.clone(),
      part_func_mat_4_rightmost_base_pairings_on_mls: neg_inf_mat.clone(),
      part_func_mat_4_base_pairings: part_func_mat.clone(),
      part_func_mat_4_base_pairings_accessible_on_el: part_func_mat.clone(),
      part_func_mat_4_base_pairings_accessible_on_mls: part_func_mat.clone(),
      part_func_mat_4_at_least_1_base_pairings_on_mls: neg_inf_mat,
      part_func_mat_4_base_unpairings_on_sa: zero_mat.clone(),
      part_func_mat_4_base_unpairings_on_el: zero_mat.clone(),
      part_func_mat_4_base_unpairings_on_mls: zero_mat,
      score_mat_4_ml_closing_basepairings: part_func_mat,
      score_mat_4_loop_aligns: vec![NEG_INFINITY; sa_len],
      twoloop_score_4d_mat: PartFunc4dMat::<T>::default(),
    }
  }
}

pub const GAP: Char = '-' as Char;
pub const CONSPROB_MAX_HAIRPIN_LOOP_LEN_CSS: usize = 2 * CONSPROB_MAX_HAIRPIN_LOOP_LEN;
pub const CONSPROB_MAX_TWOLOOP_LEN_CSS: usize = 2 * CONSPROB_MAX_TWOLOOP_LEN;

pub fn consalifold<T>(mix_bpp_mat: &ProbMat, gamma: Prob, sa: &SeqAlign<T>) -> MeaCss<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
  let sa_len = sa.cols.len();
  let mut mea_mat = vec![vec![0.; sa_len]; sa_len];
  let sa_len = T::from_usize(sa_len).unwrap();
  let gamma_plus_1 = gamma + 1.;
  for sub_sa_len in range_inclusive(T::one(), sa_len) {
    for i in range_inclusive(T::zero(), sa_len - sub_sa_len) {
      let j = i + sub_sa_len - T::one();
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
  let mut pos_pair_stack = vec![(T::zero(), sa_len - T::one())];
  while pos_pair_stack.len() > 0 {
    let pos_pair = pos_pair_stack.pop().expect("Failed to pop an element of a vector.");
    let (i, j) = pos_pair;
    if j <= i {continue;}
    let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
    let mea = mea_mat[long_i][long_j];
    if mea == mea_mat[long_i + 1][long_j] {
      pos_pair_stack.push((i + T::one(), j));
    } else if mea == mea_mat[long_i][long_j - 1] {
      pos_pair_stack.push((i, j - T::one()));
    } else if mea == mea_mat[long_i + 1][long_j - 1] + gamma_plus_1 * mix_bpp_mat[long_i][long_j] - 1. {
      pos_pair_stack.push((i + T::one(), j - T::one()));
      mea_css.bpa_pos_pairs.push(pos_pair);
    } else {
      for k in range(i, j) {
        let long_k = k.to_usize().unwrap();
        if mea == mea_mat[long_i][long_k] + mea_mat[long_k + 1][long_j] {
          pos_pair_stack.push((i, k));
          pos_pair_stack.push((k + T::one(), j));
          break;
        }
      }
    }
  }
  mea_css.ea = mea_mat[0][sa_len.to_usize().unwrap() - 1];
  mea_css
}

pub fn rnaalifold_trained<T>(sa: &SeqAlign<T>, fasta_records: &FastaRecords, feature_score_sets: &FeatureCountSets) -> SparseProbMat<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord + Sync + Send,
{
  let part_func_mats = get_css_part_func_mats(sa, fasta_records, feature_score_sets);
  get_base_pairing_prob_mat_css(sa, &part_func_mats)
}

pub fn get_css_part_func_mats<T>(sa: &SeqAlign<T>, fasta_records: &FastaRecords, feature_score_sets: &FeatureCountSets) -> CssPartFuncMats<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
  let sa_len = sa.cols.len();
  let mut css_part_func_mats = CssPartFuncMats::<T>::new(sa_len);
  for i in 0 .. sa_len {
    css_part_func_mats.score_mat_4_loop_aligns[i] = get_loop_align_score_avg(sa, i, feature_score_sets);
  }
  let short_sa_len = T::from_usize(sa_len).unwrap();
  for sub_sa_len in range_inclusive(T::one(), short_sa_len) {
    for i in range_inclusive(T::zero(), short_sa_len - sub_sa_len) {
      let j = i + sub_sa_len - T::one();
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      let pp_closing_loop = (i, j);
      let long_pp_closing_loop = (long_i, long_j);
      let loop_align_score_avg = css_part_func_mats.score_mat_4_loop_aligns[long_j];
      css_part_func_mats.part_func_mat_4_base_unpairings_on_sa[long_i][long_j] = if long_j == 0 {0.} else {css_part_func_mats.part_func_mat_4_base_unpairings_on_sa[long_i][long_j - 1]} + loop_align_score_avg;
      let num_of_nongap_nucs = sa.cols[long_j].iter().map(|&x| x != PSEUDO_BASE).len() as FeatureCount;
      css_part_func_mats.part_func_mat_4_base_unpairings_on_el[long_i][long_j] = if long_j == 0 {0.} else {css_part_func_mats.part_func_mat_4_base_unpairings_on_el[long_i][long_j - 1]} + loop_align_score_avg + num_of_nongap_nucs * feature_score_sets.external_loop_accessible_baseunpairing_count;
      css_part_func_mats.part_func_mat_4_base_unpairings_on_mls[long_i][long_j] = if long_j == 0 {0.} else {css_part_func_mats.part_func_mat_4_base_unpairings_on_mls[long_i][long_j - 1]} + loop_align_score_avg + num_of_nongap_nucs * feature_score_sets.multi_loop_accessible_baseunpairing_count;
      let mut sum = NEG_INFINITY;
      if long_pp_closing_loop.1 - long_pp_closing_loop.0 + 1 >= CONSPROB_MIN_HAIRPIN_LOOP_SPAN {
        let basepair_align_score_avg = get_basepair_align_score_avg(sa, &long_pp_closing_loop, feature_score_sets);
        if long_j - long_i - 1 <= CONSPROB_MAX_HAIRPIN_LOOP_LEN_CSS {
          let hairpin_loop_score_avg = get_hairpin_loop_score_avg(sa, fasta_records, &long_pp_closing_loop, feature_score_sets) + basepair_align_score_avg;
          logsumexp(&mut sum, hairpin_loop_score_avg);
        }
        for k in range(i + T::one(), j - T::one()) {
          let long_k = k.to_usize().unwrap();
          for l in range(k + T::one(), j) {
            let long_l = l.to_usize().unwrap();
            if long_j - long_l - 1 + long_k - long_i - 1 > CONSPROB_MAX_TWOLOOP_LEN_CSS {continue;}
            let accessible_pp = (k, l);
            let long_accessible_pp = (long_k, long_l);
            match css_part_func_mats.part_func_mat_4_base_pairings.get(&accessible_pp) {
              Some(&part_func) => {
                let twoloop_score_avg = get_twoloop_score_avg(sa, fasta_records, &long_pp_closing_loop, &long_accessible_pp, feature_score_sets) + basepair_align_score_avg + css_part_func_mats.part_func_mat_4_base_unpairings_on_sa[long_i + 1][long_k - 1] + css_part_func_mats.part_func_mat_4_base_unpairings_on_sa[long_l + 1][long_j - 1];
                css_part_func_mats.twoloop_score_4d_mat.insert((i, j, k, l), twoloop_score_avg);
                logsumexp(&mut sum, part_func + twoloop_score_avg);
              }, None => {},
            }
          }
        }
        let ml_closing_basepairing_score_avg = get_ml_closing_basepairing_score_avg(sa, fasta_records, &long_pp_closing_loop, feature_score_sets) + basepair_align_score_avg;
        css_part_func_mats.score_mat_4_ml_closing_basepairings.insert(pp_closing_loop, ml_closing_basepairing_score_avg);
        for k in long_i + 1 .. long_j {
          logsumexp(&mut sum, css_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i + 1][k - 1] + css_part_func_mats.part_func_mat_4_rightmost_base_pairings_on_mls[k][long_j - 1] + ml_closing_basepairing_score_avg);
        }
        if sum > NEG_INFINITY {
          css_part_func_mats.part_func_mat_4_base_pairings.insert(pp_closing_loop, sum);
          css_part_func_mats.part_func_mat_4_base_pairings_accessible_on_el.insert(pp_closing_loop, sum + get_el_accessible_basepairing_score_avg(sa, fasta_records, &long_pp_closing_loop, feature_score_sets));
          css_part_func_mats.part_func_mat_4_base_pairings_accessible_on_mls.insert(pp_closing_loop, sum + get_ml_accessible_basepairing_score_avg(sa, fasta_records, &long_pp_closing_loop, feature_score_sets));
        }
      }
      sum = NEG_INFINITY;
      let mut sum_2 = sum;
      for k in range_inclusive(i + T::one(), j) {
        let accessible_pp = (i, k);
        match css_part_func_mats.part_func_mat_4_base_pairings_accessible_on_el.get(&accessible_pp) {
          Some(&part_func) => {
            let long_k = k.to_usize().unwrap();
            // logsumexp(&mut sum, part_func + 2. * feature_score_sets.external_loop_accessible_baseunpairing_count * (j - k).to_f32().unwrap());
            logsumexp(&mut sum, part_func + if long_k == sa_len {0.} else {css_part_func_mats.part_func_mat_4_base_unpairings_on_el[long_k + 1][long_j]});
            // logsumexp(&mut sum_2, css_part_func_mats.part_func_mat_4_base_pairings_accessible_on_mls[&accessible_pp] + 2. * feature_score_sets.multi_loop_accessible_baseunpairing_count * (j - k).to_f32().unwrap());
            logsumexp(&mut sum_2, css_part_func_mats.part_func_mat_4_base_pairings_accessible_on_mls[&accessible_pp] + if long_k == sa_len {0.} else {css_part_func_mats.part_func_mat_4_base_unpairings_on_el[long_k + 1][long_j]});
          }, None => {},
        }
      }
      css_part_func_mats.part_func_mat_4_rightmost_base_pairings_on_el[long_i][long_j] = sum;
      css_part_func_mats.part_func_mat_4_rightmost_base_pairings_on_mls[long_i][long_j] = sum_2;
      // sum = 2. * feature_score_sets.external_loop_accessible_baseunpairing_count * sub_sa_len.to_f32().unwrap();
      sum = css_part_func_mats.part_func_mat_4_base_unpairings_on_el[long_i][long_j];
      for k in long_i .. long_j {
        let css_part_func_4_rightmost_base_pairings_on_el = css_part_func_mats.part_func_mat_4_rightmost_base_pairings_on_el[k][long_j];
        if css_part_func_4_rightmost_base_pairings_on_el == NEG_INFINITY {
          continue;
        }
        let part_func = if long_i == 0 && k == 0 {0.} else {css_part_func_mats.part_func_mat[long_i][k - 1]};
        logsumexp(&mut sum, part_func + css_part_func_4_rightmost_base_pairings_on_el);
      }
      css_part_func_mats.part_func_mat[long_i][long_j] = sum;
      sum = css_part_func_mats.part_func_mat_4_rightmost_base_pairings_on_mls[long_i][long_j];
      for k in long_i + 1 .. long_j {
        let css_part_func_4_rightmost_base_pairings_on_mls = css_part_func_mats.part_func_mat_4_rightmost_base_pairings_on_mls[k][long_j];
        // logsumexp(&mut sum, css_part_func_4_rightmost_base_pairings_on_mls + 2. * feature_score_sets.multi_loop_accessible_baseunpairing_count * (k - long_i) as FeatureCount);
        logsumexp(&mut sum, css_part_func_mats.part_func_mat_4_base_unpairings_on_mls[long_i][k - 1] + css_part_func_4_rightmost_base_pairings_on_mls);
        logsumexp(&mut sum, css_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i][k - 1] + css_part_func_4_rightmost_base_pairings_on_mls);
      }
      css_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_i][long_j] = sum;
    }
  }
  css_part_func_mats
}

pub fn get_loop_align_score_avg<T>(sa: &SeqAlign<T>, pos: usize, feature_score_sets: &FeatureCountSets) -> FeatureCount where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let mut loop_align_score_avg = 0.;
  let mut count = 0;
  let num_of_rnas = sa.cols[0].len();
  let ref col = sa.cols[pos];
  for i in 0 .. num_of_rnas {
    let base = col[i];
    for j in i + 1 .. num_of_rnas {
      let base_2 = col[j];
      if base == PSEUDO_BASE && base_2 == PSEUDO_BASE {continue;}
      loop_align_score_avg += if base != PSEUDO_BASE && base_2 != PSEUDO_BASE {
        feature_score_sets.loop_align_count_mat[base][base_2]
      } else {
        if pos == 0 {
          feature_score_sets.opening_gap_count
        } else {
          let ref prev_col = sa.cols[pos - 1];
          let prev_base_pair = (prev_col[i], prev_col[j]);
          if prev_base_pair.0 == PSEUDO_BASE && prev_base_pair.1 == PSEUDO_BASE {
            feature_score_sets.opening_gap_count
          } else if base == PSEUDO_BASE && prev_base_pair.0 != PSEUDO_BASE {
            feature_score_sets.opening_gap_count
          } else if base_2 == PSEUDO_BASE && prev_base_pair.1 != PSEUDO_BASE {
            feature_score_sets.opening_gap_count
          } else {
            feature_score_sets.extending_gap_count
          }
        }
      };
      count += 1;
    }
  }
  if count == 0 {
    NEG_INFINITY
  } else {
    loop_align_score_avg / count as FeatureCount
  }
}

fn get_base_pairing_prob_mat_css<T>(sa: &SeqAlign<T>, css_part_func_mats: &CssPartFuncMats<T>) -> SparseProbMat<T>
where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer,
{
  let sa_len = sa.cols.len();
  let css_part_func = css_part_func_mats.part_func_mat[0][sa_len - 1];
  let mut bpp_mat = SparseProbMat::<T>::default();
  let mut prob_mat_4_mls_1 = vec![vec![NEG_INFINITY; sa_len]; sa_len];
  let mut prob_mat_4_mls_2 = prob_mat_4_mls_1.clone();
  let short_sa_len = T::from_usize(sa_len).unwrap();
  for sub_sa_len in range_inclusive(T::from_usize(CONSPROB_MIN_HAIRPIN_LOOP_SPAN).unwrap(), short_sa_len).rev() {
    for i in range_inclusive(T::zero(), short_sa_len - sub_sa_len) {
      let j = i + sub_sa_len - T::one();
      let (long_i, long_j) = (i.to_usize().unwrap(), j.to_usize().unwrap());
      let mut sum_1 = NEG_INFINITY;
      let mut sum_2 = sum_1;
      for k in range(j + T::one(), short_sa_len) {
        let long_k = k.to_usize().unwrap();
        let pp_closing_loop = (i, k);
        match css_part_func_mats.part_func_mat_4_base_pairings.get(&pp_closing_loop) {
          Some(&part_func) => {
            let bpp = bpp_mat[&pp_closing_loop];
            let coefficient = bpp + css_part_func_mats.score_mat_4_ml_closing_basepairings[&pp_closing_loop] - part_func;
            logsumexp(&mut sum_1, coefficient + css_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[long_j + 1][long_k - 1]);
            // logsumexp(&mut sum_2, coefficient + 2. * feature_score_sets.multi_loop_accessible_baseunpairing_count * (k - j - T::one()).to_f32().unwrap());
            logsumexp(&mut sum_2, coefficient + css_part_func_mats.part_func_mat_4_base_unpairings_on_mls[long_j + 1][long_k - 1]);
          }, None => {},
        }
      }
      prob_mat_4_mls_1[long_i][long_j] = sum_1;
      prob_mat_4_mls_2[long_i][long_j] = sum_2;
      let accessible_pp = (i, j);
      match css_part_func_mats.part_func_mat_4_base_pairings.get(&accessible_pp) {
        Some(&part_func) => {
          let part_func_pair = (
            if accessible_pp.0 < T::one() {0.} else {css_part_func_mats.part_func_mat[0][long_i - 1]},
            if accessible_pp.1 > short_sa_len - T::from_usize(2).unwrap() {0.} else {css_part_func_mats.part_func_mat[long_j + 1][sa_len - 1]},
          );
          let mut sum = part_func_pair.0 + part_func_pair.1 + css_part_func_mats.part_func_mat_4_base_pairings_accessible_on_el[&accessible_pp] - css_part_func;
          let mut bpp_4_2l = NEG_INFINITY;
          for k in range(T::zero(), i) {
            let long_k = k.to_usize().unwrap();
            for l in range(j + T::one(), short_sa_len) {
              let long_l = l.to_usize().unwrap();
              if long_l - long_j - 1 + long_i - long_k - 1 > CONSPROB_MAX_TWOLOOP_LEN_CSS {continue;}
              let pp_closing_loop = (k, l);
              match css_part_func_mats.part_func_mat_4_base_pairings.get(&pp_closing_loop) {
                Some(&part_func_2) => {
                  logsumexp(&mut bpp_4_2l, bpp_mat[&pp_closing_loop] + part_func - part_func_2 + css_part_func_mats.twoloop_score_4d_mat[&(k, l, i, j)]);
                }, None => {},
              }
            }
          }
          if bpp_4_2l > NEG_INFINITY {
            logsumexp(&mut sum, bpp_4_2l);
          }
          let coefficient = css_part_func_mats.part_func_mat_4_base_pairings_accessible_on_mls[&accessible_pp];
          let mut bpp_4_ml = NEG_INFINITY;
          for k in 0 .. long_i {
            let css_part_func_4_at_least_1_base_pairings_on_mls = css_part_func_mats.part_func_mat_4_at_least_1_base_pairings_on_mls[k + 1][long_i - 1];
            logsumexp(&mut bpp_4_ml, coefficient + prob_mat_4_mls_2[k][long_j] + css_part_func_4_at_least_1_base_pairings_on_mls);
            let prob_4_mls = prob_mat_4_mls_1[k][long_j];
            // logsumexp(&mut bpp_4_ml, coefficient + prob_4_mls + 2. * feature_score_sets.multi_loop_accessible_baseunpairing_count * (long_i - k - 1) as FreeEnergy);
            logsumexp(&mut bpp_4_ml, coefficient + prob_4_mls + css_part_func_mats.part_func_mat_4_base_unpairings_on_mls[k + 1][long_i - 1]);
            logsumexp(&mut bpp_4_ml, coefficient + prob_4_mls + css_part_func_4_at_least_1_base_pairings_on_mls);
          }
          if bpp_4_ml > NEG_INFINITY {
            logsumexp(&mut sum, bpp_4_ml);
          }
          debug_assert!(NEG_INFINITY <= sum && sum <= 0.);
          bpp_mat.insert(accessible_pp, sum);
        }, None => {},
      }
    }
  }
  bpp_mat = bpp_mat.iter().map(|(pos_pair, &bpp)| {(*pos_pair, bpp.exp())}).collect();
  bpp_mat
}

pub fn get_hairpin_loop_score_avg<T>(sa: &SeqAlign<T>, fasta_records: &FastaRecords, pos_pair_closing_loop: &(usize, usize), feature_score_sets: &FeatureCountSets) -> FeatureCount where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let col_pair = (&sa.cols[pos_pair_closing_loop.0], &sa.cols[pos_pair_closing_loop.1]);
  let mut scores = Vec::new();
  for i in 0 .. sa.cols[0].len() {
    let ref seq = fasta_records[i].seq[..];
    let pos_pair = (sa.pos_map_sets[pos_pair_closing_loop.0][i].to_usize().unwrap(), sa.pos_map_sets[pos_pair_closing_loop.1][i].to_usize().unwrap());
    let base_pair = (col_pair.0[i], col_pair.1[i]);
    if !is_canonical(&base_pair) {continue;}
    let hairpin_loop_len = pos_pair.1 - pos_pair.0 - 1;
    if hairpin_loop_len > CONSPROB_MAX_HAIRPIN_LOOP_LEN || hairpin_loop_len < CONSPROB_MIN_HAIRPIN_LOOP_LEN {continue;}
    let hairpin_loop_score = get_consprob_hairpin_loop_score(feature_score_sets, seq, &pos_pair);
    scores.push(hairpin_loop_score);
  }
  let mut hairpin_loop_score_avg = 0.;
  let mut count = 0;
  for (i, &score) in scores.iter().enumerate() {
    for score_2 in scores[i + 1 ..].iter() {
      hairpin_loop_score_avg += score + score_2;
      count += 1;
    }
  }
  if count == 0 {
    NEG_INFINITY
  } else {
    hairpin_loop_score_avg / count as FeatureCount
  }
}

pub fn get_twoloop_score_avg<T>(
  sa: &SeqAlign<T>,
  fasta_records: &FastaRecords, 
  pos_pair_closing_loop: &(usize, usize),
  pos_pair_accessible: &(usize, usize),
  feature_score_sets: &FeatureCountSets,
) -> FeatureCount where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let col_pair = (&sa.cols[pos_pair_closing_loop.0], &sa.cols[pos_pair_closing_loop.1]);
  let col_pair_2 = (&sa.cols[pos_pair_accessible.0], &sa.cols[pos_pair_accessible.1]);
  let mut scores = Vec::new();
  for i in 0 .. sa.cols[0].len() {
    let ref seq = fasta_records[i].seq[..];
    let pos_pair = (sa.pos_map_sets[pos_pair_closing_loop.0][i].to_usize().unwrap(), sa.pos_map_sets[pos_pair_closing_loop.1][i].to_usize().unwrap());
    let base_pair = (col_pair.0[i], col_pair.1[i]);
    if !is_canonical(&base_pair) {continue;}
    let pos_pair_2 = (sa.pos_map_sets[pos_pair_accessible.0][i].to_usize().unwrap(), sa.pos_map_sets[pos_pair_accessible.1][i].to_usize().unwrap());
    let base_pair_2 = (col_pair_2.0[i], col_pair_2.1[i]);
    if !is_canonical(&base_pair_2) {continue;}
    let twoloop_len = pos_pair_2.0 - pos_pair.0 - 1 + pos_pair.1 - pos_pair_2.1 - 1;
    if twoloop_len > CONSPROB_MAX_TWOLOOP_LEN {continue;}
    let twoloop_score = get_consprob_twoloop_score(feature_score_sets, seq, &pos_pair, &pos_pair_2);
    scores.push(twoloop_score);
  }
  let mut twoloop_score_avg = 0.;
  let mut count = 0;
  for (i, &score) in scores.iter().enumerate() {
    for score_2 in scores[i + 1 ..].iter() {
      twoloop_score_avg += score + score_2;
      count += 1;
    }
  }
  if count == 0 {
    NEG_INFINITY
  } else {
    twoloop_score_avg / count as FeatureCount
  }
}

pub fn get_ml_closing_basepairing_score_avg<T>(sa: &SeqAlign<T>, fasta_records: &FastaRecords, pos_pair_closing_loop: &(usize, usize), feature_score_sets: &FeatureCountSets) -> FeatureCount where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let col_pair = (&sa.cols[pos_pair_closing_loop.0], &sa.cols[pos_pair_closing_loop.1]);
  let mut scores = Vec::new();
  for i in 0 .. sa.cols[0].len() {
    let ref seq = fasta_records[i].seq[..];
    let pos_pair = (sa.pos_map_sets[pos_pair_closing_loop.0][i].to_usize().unwrap(), sa.pos_map_sets[pos_pair_closing_loop.1][i].to_usize().unwrap());
    let base_pair = (col_pair.0[i], col_pair.1[i]);
    if !is_canonical(&base_pair) {continue;}
    let multi_loop_closing_basepairing_score = get_consprob_multi_loop_closing_basepairing_score(feature_score_sets, seq, &pos_pair);
    scores.push(multi_loop_closing_basepairing_score);
  }
  let mut ml_closing_basepairing_score_avg = 0.;
  let mut count = 0;
  for (i, &score) in scores.iter().enumerate() {
    for score_2 in scores[i + 1 ..].iter() {
      ml_closing_basepairing_score_avg += score + score_2;
      count += 1;
    }
  }
  if count == 0 {
    NEG_INFINITY
  } else {
    ml_closing_basepairing_score_avg / count as FeatureCount
  }
}

pub fn get_el_accessible_basepairing_score_avg<T>(sa: &SeqAlign<T>, fasta_records: &FastaRecords, pos_pair_accessible: &(usize, usize), feature_score_sets: &FeatureCountSets) -> FeatureCount where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let col_pair = (&sa.cols[pos_pair_accessible.0], &sa.cols[pos_pair_accessible.1]);
  let mut scores = Vec::new();
  for i in 0 .. sa.cols[0].len() {
    let ref seq = fasta_records[i].seq[..];
    let pos_pair = (sa.pos_map_sets[pos_pair_accessible.0][i].to_usize().unwrap(), sa.pos_map_sets[pos_pair_accessible.1][i].to_usize().unwrap());
    let base_pair = (col_pair.0[i], col_pair.1[i]);
    if !is_canonical(&base_pair) {continue;}
    let external_loop_accessible_basepairing_score = get_consprob_external_loop_accessible_basepairing_score(feature_score_sets, seq, &pos_pair);
    scores.push(external_loop_accessible_basepairing_score);
  }
  let mut el_accessible_basepairing_score_avg = 0.;
  let mut count = 0;
  for (i, &score) in scores.iter().enumerate() {
    for score_2 in scores[i + 1 ..].iter() {
      el_accessible_basepairing_score_avg += score + score_2;
      count += 1;
    }
  }
  if count == 0 {
    NEG_INFINITY
  } else {
    el_accessible_basepairing_score_avg / count as FeatureCount
  }
}

pub fn get_ml_accessible_basepairing_score_avg<T>(sa: &SeqAlign<T>, fasta_records: &FastaRecords, pos_pair_accessible: &(usize, usize), feature_score_sets: &FeatureCountSets) -> FeatureCount where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let col_pair = (&sa.cols[pos_pair_accessible.0], &sa.cols[pos_pair_accessible.1]);
  let mut scores = Vec::new();
  for i in 0 .. sa.cols[0].len() {
    let ref seq = fasta_records[i].seq[..];
    let pos_pair = (sa.pos_map_sets[pos_pair_accessible.0][i].to_usize().unwrap(), sa.pos_map_sets[pos_pair_accessible.1][i].to_usize().unwrap());
    let base_pair = (col_pair.0[i], col_pair.1[i]);
    if !is_canonical(&base_pair) {continue;}
    let multi_loop_accessible_basepairing_score = get_consprob_multi_loop_accessible_basepairing_score(feature_score_sets, seq, &pos_pair);
    scores.push(multi_loop_accessible_basepairing_score);
  }
  let mut ml_accessible_basepairing_score_avg = 0.;
  let mut count = 0;
  for (i, &score) in scores.iter().enumerate() {
    for score_2 in scores[i + 1 ..].iter() {
      ml_accessible_basepairing_score_avg += score + score_2;
      count += 1;
    }
  }
  if count == 0 {
    NEG_INFINITY
  } else {
    ml_accessible_basepairing_score_avg / count as FeatureCount
  }
}

pub fn get_basepair_align_score_avg<T>(sa: &SeqAlign<T>, pos_pair: &(usize, usize), feature_score_sets: &FeatureCountSets) -> FeatureCount where
  T: Unsigned + PrimInt + Hash + FromPrimitive + Integer + Ord,
{
  let mut basepair_align_score_avg = 0.;
  let mut count = 0;
  let num_of_rnas = sa.cols[0].len();
  let col_pair = (&sa.cols[pos_pair.0], &sa.cols[pos_pair.1]);
  for i in 0 .. num_of_rnas {
    let base_pair = (col_pair.0[i], col_pair.1[i]);
    if base_pair.0 == PSEUDO_BASE || base_pair.1 == PSEUDO_BASE {continue;}
    for j in i + 1 .. num_of_rnas {
      let base_pair_2 = (col_pair.0[j], col_pair.1[j]);
      if base_pair_2.0 == PSEUDO_BASE || base_pair_2.1 == PSEUDO_BASE {continue;}
      let basepair_align_score = feature_score_sets.basepair_align_count_mat[base_pair.0][base_pair.1][base_pair_2.0][base_pair_2.1];
      basepair_align_score_avg += basepair_align_score;
      count += 1;
    }
  }
  if count == 0 {
    NEG_INFINITY
  } else {
    basepair_align_score_avg / count as FeatureCount
  }
}

pub fn revert_char(c: Base) -> u8 {
  match c {
    A => BIG_A,
    C => BIG_C,
    G => BIG_G,
    U => BIG_U,
    PSEUDO_BASE => {GAP},
    _ => {assert!(false); GAP},
  }
}
