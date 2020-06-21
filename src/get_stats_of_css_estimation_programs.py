#! /usr/bin/env python

import utils
from Bio import SeqIO
import numpy
import seaborn
from matplotlib import pyplot
import os
from sklearn.metrics import roc_curve
from sklearn.metrics import auc
import math

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  short_mafft_ginsi_plus_consalifold_ppvs = []
  short_mafft_ginsi_plus_consalifold_senss = []
  short_mafft_ginsi_plus_consalifold_fprs = []
  short_mafft_xinsi_plus_consalifold_ppvs = []
  short_mafft_xinsi_plus_consalifold_senss = []
  short_mafft_xinsi_plus_consalifold_fprs = []
  short_ref_sa_plus_consalifold_ppvs = []
  short_ref_sa_plus_consalifold_senss = []
  short_ref_sa_plus_consalifold_fprs = []
  short_mafft_ginsi_plus_centroidalifold_ppvs = []
  short_mafft_ginsi_plus_centroidalifold_senss = []
  short_mafft_ginsi_plus_centroidalifold_fprs = []
  short_mafft_xinsi_plus_centroidalifold_ppvs = []
  short_mafft_xinsi_plus_centroidalifold_senss = []
  short_mafft_xinsi_plus_centroidalifold_fprs = []
  short_ref_sa_plus_centroidalifold_ppvs = []
  short_ref_sa_plus_centroidalifold_senss = []
  short_ref_sa_plus_centroidalifold_fprs = []
  short_mafft_ginsi_plus_rnaalifold_ppv = short_mafft_ginsi_plus_rnaalifold_sens = short_mafft_ginsi_plus_rnaalifold_fpr = 0.
  short_mafft_xinsi_plus_rnaalifold_ppv = short_mafft_xinsi_plus_rnaalifold_sens = short_mafft_xinsi_plus_rnaalifold_fpr = 0.
  short_ref_sa_plus_rnaalifold_ppv = short_ref_sa_plus_rnaalifold_sens = short_ref_sa_plus_rnaalifold_fpr = 0.
  long_mafft_ginsi_plus_consalifold_ppvs = []
  long_mafft_ginsi_plus_consalifold_senss = []
  long_mafft_ginsi_plus_consalifold_fprs = []
  long_mafft_xinsi_plus_consalifold_ppvs = []
  long_mafft_xinsi_plus_consalifold_senss = []
  long_mafft_xinsi_plus_consalifold_fprs = []
  long_ref_sa_plus_consalifold_ppvs = []
  long_ref_sa_plus_consalifold_senss = []
  long_ref_sa_plus_consalifold_fprs = []
  long_mafft_ginsi_plus_centroidalifold_ppvs = []
  long_mafft_ginsi_plus_centroidalifold_senss = []
  long_mafft_ginsi_plus_centroidalifold_fprs = []
  long_mafft_xinsi_plus_centroidalifold_ppvs = []
  long_mafft_xinsi_plus_centroidalifold_senss = []
  long_mafft_xinsi_plus_centroidalifold_fprs = []
  long_ref_sa_plus_centroidalifold_ppvs = []
  long_ref_sa_plus_centroidalifold_senss = []
  long_ref_sa_plus_centroidalifold_fprs = []
  long_mafft_ginsi_plus_rnaalifold_ppv = long_mafft_ginsi_plus_rnaalifold_sens = long_mafft_ginsi_plus_rnaalifold_fpr = 0.
  long_mafft_xinsi_plus_rnaalifold_ppv = long_mafft_xinsi_plus_rnaalifold_sens = long_mafft_xinsi_plus_rnaalifold_fpr = 0.
  long_ref_sa_plus_rnaalifold_ppv = long_ref_sa_plus_rnaalifold_sens = long_ref_sa_plus_rnaalifold_fpr = 0.
  gammas = [2. ** i for i in range(-7, 11)]
  for data_set in ["short", "long"]:
    # rna_fam_dir_path = asset_dir_path + "/compiled_rna_fams_" + data_set
    rna_fam_dir_path = asset_dir_path + "/compiled_rna_fams_" + data_set + "_4_micro_bench"
    # ref_sa_dir_path = asset_dir_path + "/ref_sas_" + data_set
    ref_sa_dir_path = asset_dir_path + "/ref_sas_" + data_set + "_4_micro_bench"
    mafft_ginsi_plus_consalifold_css_dir_path = asset_dir_path + "/mafft_ginsi_plus_consalifold_" + data_set
    mafft_xinsi_plus_consalifold_css_dir_path = asset_dir_path + "/mafft_xinsi_plus_consalifold_" + data_set
    ref_sa_plus_consalifold_css_dir_path = asset_dir_path + "/ref_sa_plus_consalifold_" + data_set
    mafft_ginsi_plus_centroidalifold_css_dir_path = asset_dir_path + "/mafft_ginsi_plus_centroidalifold_" + data_set
    mafft_xinsi_plus_centroidalifold_css_dir_path = asset_dir_path + "/mafft_xinsi_plus_centroidalifold_" + data_set
    ref_sa_plus_centroidalifold_css_dir_path = asset_dir_path + "/ref_sa_plus_centroidalifold_" + data_set
    mafft_ginsi_plus_rnaalifold_css_dir_path = asset_dir_path + "/mafft_ginsi_plus_rnaalifold_" + data_set
    mafft_xinsi_plus_rnaalifold_css_dir_path = asset_dir_path + "/mafft_xinsi_plus_rnaalifold_" + data_set
    ref_sa_plus_rnaalifold_css_dir_path = asset_dir_path + "/ref_sa_plus_rnaalifold_" + data_set
    for gamma in gammas:
      gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
      mafft_ginsi_plus_consalifold_tp = mafft_ginsi_plus_consalifold_tn = mafft_ginsi_plus_consalifold_fp = mafft_ginsi_plus_consalifold_fn = 0.
      mafft_xinsi_plus_consalifold_tp = mafft_xinsi_plus_consalifold_tn = mafft_xinsi_plus_consalifold_fp = mafft_xinsi_plus_consalifold_fn = 0.
      ref_sa_plus_consalifold_tp = ref_sa_plus_consalifold_tn = ref_sa_plus_consalifold_fp = ref_sa_plus_consalifold_fn = 0.
      mafft_ginsi_plus_centroidalifold_tp = mafft_ginsi_plus_centroidalifold_tn = mafft_ginsi_plus_centroidalifold_fp = mafft_ginsi_plus_centroidalifold_fn = 0.
      mafft_xinsi_plus_centroidalifold_tp = mafft_xinsi_plus_centroidalifold_tn = mafft_xinsi_plus_centroidalifold_fp = mafft_xinsi_plus_centroidalifold_fn = 0.
      ref_sa_plus_centroidalifold_tp = ref_sa_plus_centroidalifold_tn = ref_sa_plus_centroidalifold_fp = ref_sa_plus_centroidalifold_fn = 0.
      mafft_ginsi_plus_rnaalifold_tp = mafft_ginsi_plus_rnaalifold_tn = mafft_ginsi_plus_rnaalifold_fp = mafft_ginsi_plus_rnaalifold_fn = 0.
      mafft_xinsi_plus_rnaalifold_tp = mafft_xinsi_plus_rnaalifold_tn = mafft_xinsi_plus_rnaalifold_fp = mafft_xinsi_plus_rnaalifold_fn = 0.
      ref_sa_plus_rnaalifold_tp = ref_sa_plus_rnaalifold_tn = ref_sa_plus_rnaalifold_fp = ref_sa_plus_rnaalifold_fn = 0.
      for rna_fam_file in os.listdir(rna_fam_dir_path):
        if not rna_fam_file.endswith(".fa"):
          continue
        rna_seq_file_path = os.path.join(rna_fam_dir_path, rna_fam_file)
        rna_seq_lens = [len(rna_seq.seq) for rna_seq in SeqIO.parse(rna_seq_file_path, "fasta")]
        num_of_rnas = len(rna_seq_lens)
        (rna_fam_name, extension) = os.path.splitext(rna_fam_file)
        ref_css_file_path = os.path.join(ref_sa_dir_path, rna_fam_name + ".sth")
        ref_css_and_flat_css = utils.get_css_and_flat_css(ref_css_file_path)
        mafft_ginsi_plus_consalifold_estimated_css_dir_path = os.path.join(mafft_ginsi_plus_consalifold_css_dir_path, rna_fam_name)
        if not os.path.isdir(mafft_ginsi_plus_consalifold_estimated_css_dir_path):
          continue
        mafft_xinsi_plus_consalifold_estimated_css_dir_path = os.path.join(mafft_xinsi_plus_consalifold_css_dir_path, rna_fam_name)
        if not os.path.isdir(mafft_xinsi_plus_consalifold_estimated_css_dir_path):
          continue
        ref_sa_plus_consalifold_estimated_css_dir_path = os.path.join(ref_sa_plus_consalifold_css_dir_path, rna_fam_name)
        if not os.path.isdir(ref_sa_plus_consalifold_estimated_css_dir_path):
          continue
        mafft_ginsi_plus_centroidalifold_estimated_css_dir_path = os.path.join(mafft_ginsi_plus_centroidalifold_css_dir_path, rna_fam_name)
        if not os.path.isdir(mafft_ginsi_plus_centroidalifold_estimated_css_dir_path):
          continue
        mafft_xinsi_plus_centroidalifold_estimated_css_dir_path = os.path.join(mafft_xinsi_plus_centroidalifold_css_dir_path, rna_fam_name)
        if not os.path.isdir(mafft_xinsi_plus_centroidalifold_estimated_css_dir_path):
          continue
        ref_sa_plus_centroidalifold_estimated_css_dir_path = os.path.join(ref_sa_plus_centroidalifold_css_dir_path, rna_fam_name)
        if not os.path.isdir(ref_sa_plus_centroidalifold_estimated_css_dir_path):
          continue
        mafft_ginsi_plus_consalifold_estimated_css_file_path = os.path.join(mafft_ginsi_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
        estimated_css_and_flat_css = utils.get_css_and_flat_css(mafft_ginsi_plus_consalifold_estimated_css_file_path)
        tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
        mafft_ginsi_plus_consalifold_tp += tp
        mafft_ginsi_plus_consalifold_tn += tn
        mafft_ginsi_plus_consalifold_fp += fp
        mafft_ginsi_plus_consalifold_fn += fn
        mafft_xinsi_plus_consalifold_estimated_css_file_path = os.path.join(mafft_xinsi_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
        estimated_css_and_flat_css = utils.get_css_and_flat_css(mafft_xinsi_plus_consalifold_estimated_css_file_path)
        tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
        mafft_xinsi_plus_consalifold_tp += tp
        mafft_xinsi_plus_consalifold_tn += tn
        mafft_xinsi_plus_consalifold_fp += fp
        mafft_xinsi_plus_consalifold_fn += fn
        ref_sa_plus_consalifold_estimated_css_file_path = os.path.join(ref_sa_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
        estimated_css_and_flat_css = utils.get_css_and_flat_css(ref_sa_plus_consalifold_estimated_css_file_path)
        tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
        ref_sa_plus_consalifold_tp += tp
        ref_sa_plus_consalifold_tn += tn
        ref_sa_plus_consalifold_fp += fp
        ref_sa_plus_consalifold_fn += fn
        mafft_ginsi_plus_centroidalifold_estimated_css_file_path = os.path.join(mafft_ginsi_plus_centroidalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
        estimated_css_and_flat_css = utils.get_css_and_flat_css(mafft_ginsi_plus_centroidalifold_estimated_css_file_path)
        tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
        mafft_ginsi_plus_centroidalifold_tp += tp
        mafft_ginsi_plus_centroidalifold_tn += tn
        mafft_ginsi_plus_centroidalifold_fp += fp
        mafft_ginsi_plus_centroidalifold_fn += fn
        mafft_xinsi_plus_centroidalifold_estimated_css_file_path = os.path.join(mafft_xinsi_plus_centroidalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
        estimated_css_and_flat_css = utils.get_css_and_flat_css(mafft_xinsi_plus_centroidalifold_estimated_css_file_path)
        tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
        mafft_xinsi_plus_centroidalifold_tp += tp
        mafft_xinsi_plus_centroidalifold_tn += tn
        mafft_xinsi_plus_centroidalifold_fp += fp
        mafft_xinsi_plus_centroidalifold_fn += fn
        ref_sa_plus_centroidalifold_estimated_css_file_path = os.path.join(ref_sa_plus_centroidalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
        estimated_css_and_flat_css = utils.get_css_and_flat_css(ref_sa_plus_centroidalifold_estimated_css_file_path)
        tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
        ref_sa_plus_centroidalifold_tp += tp
        ref_sa_plus_centroidalifold_tn += tn
        ref_sa_plus_centroidalifold_fp += fp
        ref_sa_plus_centroidalifold_fn += fn
        if gamma == 4.:
          mafft_ginsi_plus_rnaalifold_estimated_css_file_path = os.path.join(mafft_ginsi_plus_rnaalifold_css_dir_path, rna_fam_name + ".sth")
          estimated_css_and_flat_css = utils.get_css_and_flat_css(mafft_ginsi_plus_rnaalifold_estimated_css_file_path)
          tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
          mafft_ginsi_plus_rnaalifold_tp += tp
          mafft_ginsi_plus_rnaalifold_tn += tn
          mafft_ginsi_plus_rnaalifold_fp += fp
          mafft_ginsi_plus_rnaalifold_fn += fn
          mafft_xinsi_plus_rnaalifold_estimated_css_file_path = os.path.join(mafft_xinsi_plus_rnaalifold_css_dir_path, rna_fam_name + ".sth")
          estimated_css_and_flat_css = utils.get_css_and_flat_css(mafft_xinsi_plus_rnaalifold_estimated_css_file_path)
          tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
          mafft_xinsi_plus_rnaalifold_tp += tp
          mafft_xinsi_plus_rnaalifold_tn += tn
          mafft_xinsi_plus_rnaalifold_fp += fp
          mafft_xinsi_plus_rnaalifold_fn += fn
          ref_sa_plus_rnaalifold_estimated_css_file_path = os.path.join(ref_sa_plus_rnaalifold_css_dir_path, rna_fam_name + ".sth")
          estimated_css_and_flat_css = utils.get_css_and_flat_css(ref_sa_plus_rnaalifold_estimated_css_file_path)
          tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
          ref_sa_plus_rnaalifold_tp += tp
          ref_sa_plus_rnaalifold_tn += tn
          ref_sa_plus_rnaalifold_fp += fp
          ref_sa_plus_rnaalifold_fn += fn
      ppv = mafft_ginsi_plus_consalifold_tp / (mafft_ginsi_plus_consalifold_tp + mafft_ginsi_plus_consalifold_fp)
      sens = mafft_ginsi_plus_consalifold_tp / (mafft_ginsi_plus_consalifold_tp + mafft_ginsi_plus_consalifold_fn)
      fpr = mafft_ginsi_plus_consalifold_fp / (mafft_ginsi_plus_consalifold_tn + mafft_ginsi_plus_consalifold_fp)
      if data_set == "short":
        short_mafft_ginsi_plus_consalifold_ppvs.insert(0, ppv)
        short_mafft_ginsi_plus_consalifold_senss.insert(0, sens)
        short_mafft_ginsi_plus_consalifold_fprs.insert(0, fpr)
      else:
        long_mafft_ginsi_plus_consalifold_ppvs.insert(0, ppv)
        long_mafft_ginsi_plus_consalifold_senss.insert(0, sens)
        long_mafft_ginsi_plus_consalifold_fprs.insert(0, fpr)
      ppv = mafft_xinsi_plus_consalifold_tp / (mafft_xinsi_plus_consalifold_tp + mafft_xinsi_plus_consalifold_fp)
      sens = mafft_xinsi_plus_consalifold_tp / (mafft_xinsi_plus_consalifold_tp + mafft_xinsi_plus_consalifold_fn)
      fpr = mafft_xinsi_plus_consalifold_fp / (mafft_xinsi_plus_consalifold_tn + mafft_xinsi_plus_consalifold_fp)
      if data_set == "short":
        short_mafft_xinsi_plus_consalifold_ppvs.insert(0, ppv)
        short_mafft_xinsi_plus_consalifold_senss.insert(0, sens)
        short_mafft_xinsi_plus_consalifold_fprs.insert(0, fpr)
      else:
        long_mafft_xinsi_plus_consalifold_ppvs.insert(0, ppv)
        long_mafft_xinsi_plus_consalifold_senss.insert(0, sens)
        long_mafft_xinsi_plus_consalifold_fprs.insert(0, fpr)
      ppv = ref_sa_plus_consalifold_tp / (ref_sa_plus_consalifold_tp + ref_sa_plus_consalifold_fp)
      sens = ref_sa_plus_consalifold_tp / (ref_sa_plus_consalifold_tp + ref_sa_plus_consalifold_fn)
      fpr = ref_sa_plus_consalifold_fp / (ref_sa_plus_consalifold_tn + ref_sa_plus_consalifold_fp)
      if data_set == "short":
        short_ref_sa_plus_consalifold_ppvs.insert(0, ppv)
        short_ref_sa_plus_consalifold_senss.insert(0, sens)
        short_ref_sa_plus_consalifold_fprs.insert(0, fpr)
      else:
        long_ref_sa_plus_consalifold_ppvs.insert(0, ppv)
        long_ref_sa_plus_consalifold_senss.insert(0, sens)
        long_ref_sa_plus_consalifold_fprs.insert(0, fpr)
      ppv = mafft_ginsi_plus_centroidalifold_tp / (mafft_ginsi_plus_centroidalifold_tp + mafft_ginsi_plus_centroidalifold_fp)
      sens = mafft_ginsi_plus_centroidalifold_tp / (mafft_ginsi_plus_centroidalifold_tp + mafft_ginsi_plus_centroidalifold_fn)
      fpr = mafft_ginsi_plus_centroidalifold_fp / (mafft_ginsi_plus_centroidalifold_tn + mafft_ginsi_plus_centroidalifold_fp)
      if data_set == "short":
        short_mafft_ginsi_plus_centroidalifold_ppvs.insert(0, ppv)
        short_mafft_ginsi_plus_centroidalifold_senss.insert(0, sens)
        short_mafft_ginsi_plus_centroidalifold_fprs.insert(0, fpr)
      else:
        long_mafft_ginsi_plus_centroidalifold_ppvs.insert(0, ppv)
        long_mafft_ginsi_plus_centroidalifold_senss.insert(0, sens)
        long_mafft_ginsi_plus_centroidalifold_fprs.insert(0, fpr)
      ppv = mafft_xinsi_plus_centroidalifold_tp / (mafft_xinsi_plus_centroidalifold_tp + mafft_xinsi_plus_centroidalifold_fp)
      sens = mafft_xinsi_plus_centroidalifold_tp / (mafft_xinsi_plus_centroidalifold_tp + mafft_xinsi_plus_centroidalifold_fn)
      fpr = mafft_xinsi_plus_centroidalifold_fp / (mafft_xinsi_plus_centroidalifold_tn + mafft_xinsi_plus_centroidalifold_fp)
      if data_set == "short":
        short_mafft_xinsi_plus_centroidalifold_ppvs.insert(0, ppv)
        short_mafft_xinsi_plus_centroidalifold_senss.insert(0, sens)
        short_mafft_xinsi_plus_centroidalifold_fprs.insert(0, fpr)
      else:
        long_mafft_xinsi_plus_centroidalifold_ppvs.insert(0, ppv)
        long_mafft_xinsi_plus_centroidalifold_senss.insert(0, sens)
        long_mafft_xinsi_plus_centroidalifold_fprs.insert(0, fpr)
      ppv = ref_sa_plus_centroidalifold_tp / (ref_sa_plus_centroidalifold_tp + ref_sa_plus_centroidalifold_fp)
      sens = ref_sa_plus_centroidalifold_tp / (ref_sa_plus_centroidalifold_tp + ref_sa_plus_centroidalifold_fn)
      fpr = ref_sa_plus_centroidalifold_fp / (ref_sa_plus_centroidalifold_tn + ref_sa_plus_centroidalifold_fp)
      if data_set == "short":
        short_ref_sa_plus_centroidalifold_ppvs.insert(0, ppv)
        short_ref_sa_plus_centroidalifold_senss.insert(0, sens)
        short_ref_sa_plus_centroidalifold_fprs.insert(0, fpr)
      else:
        long_ref_sa_plus_centroidalifold_ppvs.insert(0, ppv)
        long_ref_sa_plus_centroidalifold_senss.insert(0, sens)
        long_ref_sa_plus_centroidalifold_fprs.insert(0, fpr)
      if gamma == 4.:
        ppv = mafft_ginsi_plus_rnaalifold_tp / (mafft_ginsi_plus_rnaalifold_tp + mafft_ginsi_plus_rnaalifold_fp)
        sens = mafft_ginsi_plus_rnaalifold_tp / (mafft_ginsi_plus_rnaalifold_tp + mafft_ginsi_plus_rnaalifold_fn)
        fpr = mafft_ginsi_plus_rnaalifold_fp / (mafft_ginsi_plus_rnaalifold_tn + mafft_ginsi_plus_rnaalifold_fp)
        if data_set == "short":
          short_mafft_ginsi_plus_rnaalifold_ppv = ppv
          short_mafft_ginsi_plus_rnaalifold_sens = sens
          short_mafft_ginsi_plus_rnaalifold_fpr = fpr
        else:
          long_mafft_ginsi_plus_rnaalifold_ppv = ppv
          long_mafft_ginsi_plus_rnaalifold_sens = sens
          long_mafft_ginsi_plus_rnaalifold_fpr = fpr
        ppv = mafft_xinsi_plus_rnaalifold_tp / (mafft_xinsi_plus_rnaalifold_tp + mafft_xinsi_plus_rnaalifold_fp)
        sens = mafft_xinsi_plus_rnaalifold_tp / (mafft_xinsi_plus_rnaalifold_tp + mafft_xinsi_plus_rnaalifold_fn)
        fpr = mafft_xinsi_plus_rnaalifold_fp / (mafft_xinsi_plus_rnaalifold_tn + mafft_xinsi_plus_rnaalifold_fp)
        if data_set == "short":
          short_mafft_xinsi_plus_rnaalifold_ppv = ppv
          short_mafft_xinsi_plus_rnaalifold_sens = sens
          short_mafft_xinsi_plus_rnaalifold_fpr = fpr
        else:
          long_mafft_xinsi_plus_rnaalifold_ppv = ppv
          long_mafft_xinsi_plus_rnaalifold_sens = sens
          long_mafft_xinsi_plus_rnaalifold_fpr = fpr
        ppv = ref_sa_plus_rnaalifold_tp / (ref_sa_plus_rnaalifold_tp + ref_sa_plus_rnaalifold_fp)
        sens = ref_sa_plus_rnaalifold_tp / (ref_sa_plus_rnaalifold_tp + ref_sa_plus_rnaalifold_fn)
        fpr = ref_sa_plus_rnaalifold_fp / (ref_sa_plus_rnaalifold_tn + ref_sa_plus_rnaalifold_fp)
        if data_set == "short":
          short_ref_sa_plus_rnaalifold_ppv = ppv
          short_ref_sa_plus_rnaalifold_sens = sens
          short_ref_sa_plus_rnaalifold_fpr = fpr
        else:
          long_ref_sa_plus_rnaalifold_ppv = ppv
          long_ref_sa_plus_rnaalifold_sens = sens
          long_ref_sa_plus_rnaalifold_fpr = fpr
  short_mafft_ginsi_plus_consalifold_ppvs = numpy.array(short_mafft_ginsi_plus_consalifold_ppvs)
  short_mafft_ginsi_plus_consalifold_senss = numpy.array(short_mafft_ginsi_plus_consalifold_senss)
  short_mafft_ginsi_plus_consalifold_fprs = numpy.array(short_mafft_ginsi_plus_consalifold_fprs)
  short_mafft_xinsi_plus_consalifold_ppvs = numpy.array(short_mafft_xinsi_plus_consalifold_ppvs)
  short_mafft_xinsi_plus_consalifold_senss = numpy.array(short_mafft_xinsi_plus_consalifold_senss)
  short_mafft_xinsi_plus_consalifold_fprs = numpy.array(short_mafft_xinsi_plus_consalifold_fprs)
  short_ref_sa_plus_consalifold_ppvs = numpy.array(short_ref_sa_plus_consalifold_ppvs)
  short_ref_sa_plus_consalifold_senss = numpy.array(short_ref_sa_plus_consalifold_senss)
  short_ref_sa_plus_consalifold_fprs = numpy.array(short_ref_sa_plus_consalifold_fprs)
  short_mafft_ginsi_plus_centroidalifold_ppvs = numpy.array(short_mafft_ginsi_plus_centroidalifold_ppvs)
  short_mafft_ginsi_plus_centroidalifold_senss = numpy.array(short_mafft_ginsi_plus_centroidalifold_senss)
  short_mafft_ginsi_plus_centroidalifold_fprs = numpy.array(short_mafft_ginsi_plus_centroidalifold_fprs)
  short_mafft_xinsi_plus_centroidalifold_ppvs = numpy.array(short_mafft_xinsi_plus_centroidalifold_ppvs)
  short_mafft_xinsi_plus_centroidalifold_senss = numpy.array(short_mafft_xinsi_plus_centroidalifold_senss)
  short_mafft_xinsi_plus_centroidalifold_fprs = numpy.array(short_mafft_xinsi_plus_centroidalifold_fprs)
  short_ref_sa_plus_centroidalifold_ppvs = numpy.array(short_ref_sa_plus_centroidalifold_ppvs)
  short_ref_sa_plus_centroidalifold_senss = numpy.array(short_ref_sa_plus_centroidalifold_senss)
  short_ref_sa_plus_centroidalifold_fprs = numpy.array(short_ref_sa_plus_centroidalifold_fprs)
  long_mafft_ginsi_plus_consalifold_ppvs = numpy.array(long_mafft_ginsi_plus_consalifold_ppvs)
  long_mafft_ginsi_plus_consalifold_senss = numpy.array(long_mafft_ginsi_plus_consalifold_senss)
  long_mafft_ginsi_plus_consalifold_fprs = numpy.array(long_mafft_ginsi_plus_consalifold_fprs)
  long_mafft_xinsi_plus_consalifold_ppvs = numpy.array(long_mafft_xinsi_plus_consalifold_ppvs)
  long_mafft_xinsi_plus_consalifold_senss = numpy.array(long_mafft_xinsi_plus_consalifold_senss)
  long_mafft_xinsi_plus_consalifold_fprs = numpy.array(long_mafft_xinsi_plus_consalifold_fprs)
  long_ref_sa_plus_consalifold_ppvs = numpy.array(long_ref_sa_plus_consalifold_ppvs)
  long_ref_sa_plus_consalifold_senss = numpy.array(long_ref_sa_plus_consalifold_senss)
  long_ref_sa_plus_consalifold_fprs = numpy.array(long_ref_sa_plus_consalifold_fprs)
  long_mafft_ginsi_plus_centroidalifold_ppvs = numpy.array(long_mafft_ginsi_plus_centroidalifold_ppvs)
  long_mafft_ginsi_plus_centroidalifold_senss = numpy.array(long_mafft_ginsi_plus_centroidalifold_senss)
  long_mafft_ginsi_plus_centroidalifold_fprs = numpy.array(long_mafft_ginsi_plus_centroidalifold_fprs)
  long_mafft_xinsi_plus_centroidalifold_ppvs = numpy.array(long_mafft_xinsi_plus_centroidalifold_ppvs)
  long_mafft_xinsi_plus_centroidalifold_senss = numpy.array(long_mafft_xinsi_plus_centroidalifold_senss)
  long_mafft_xinsi_plus_centroidalifold_fprs = numpy.array(long_mafft_xinsi_plus_centroidalifold_fprs)
  long_ref_sa_plus_centroidalifold_ppvs = numpy.array(long_ref_sa_plus_centroidalifold_ppvs)
  long_ref_sa_plus_centroidalifold_senss = numpy.array(long_ref_sa_plus_centroidalifold_senss)
  long_ref_sa_plus_centroidalifold_fprs = numpy.array(long_ref_sa_plus_centroidalifold_fprs)
  line_1, = pyplot.plot(short_mafft_ginsi_plus_consalifold_ppvs, short_mafft_ginsi_plus_consalifold_senss, label = "MAFFT G-INS-i + ConsAlifold (Short)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(short_mafft_xinsi_plus_consalifold_ppvs, short_mafft_xinsi_plus_consalifold_senss, label = "MAFFT X-INS-i + ConsAlifold (Short)", marker = "o", linestyle = "--")
  line_3, = pyplot.plot(short_ref_sa_plus_consalifold_ppvs, short_ref_sa_plus_consalifold_senss, label = "Reference + ConsAlifold (Short)", marker = "o", linestyle = ":")
  line_4, = pyplot.plot(short_mafft_ginsi_plus_centroidalifold_ppvs, short_mafft_ginsi_plus_centroidalifold_senss, label = "MAFFT G-INS-i + CentroidAlifold (Short)", marker = "s", linestyle = "-")
  line_5, = pyplot.plot(short_mafft_xinsi_plus_centroidalifold_ppvs, short_mafft_xinsi_plus_centroidalifold_senss, label = "MAFFT X-INS-i + CentroidAlifold (Short)", marker = "s", linestyle = "--")
  line_6, = pyplot.plot(short_ref_sa_plus_centroidalifold_ppvs, short_ref_sa_plus_centroidalifold_senss, label = "Reference + CentroidAlifold (Short)", marker = "s", linestyle = ":")
  line_7, = pyplot.plot(short_mafft_ginsi_plus_rnaalifold_ppv, short_mafft_ginsi_plus_rnaalifold_sens, label = "MAFFT G-INS-i + RNAalifold (Short)", marker = "v", linestyle = "-")
  line_8, = pyplot.plot(short_mafft_xinsi_plus_rnaalifold_ppv, short_mafft_xinsi_plus_rnaalifold_sens, label = "MAFFT X-INS-i + RNAalifold (Short)", marker = "v", linestyle = "--")
  line_9, = pyplot.plot(short_ref_sa_plus_rnaalifold_ppv, short_ref_sa_plus_rnaalifold_sens, label = "Reference + RNAalifold (Short)", marker = "v", linestyle = ":")
  line_10, = pyplot.plot(long_mafft_ginsi_plus_consalifold_ppvs, long_mafft_ginsi_plus_consalifold_senss, label = "MAFFT G-INS-i + ConsAlifold (Long)", marker = "^", linestyle = "-")
  line_11, = pyplot.plot(long_mafft_xinsi_plus_consalifold_ppvs, long_mafft_xinsi_plus_consalifold_senss, label = "MAFFT X-INS-i + ConsAlifold (Long)", marker = "^", linestyle = "--")
  line_12, = pyplot.plot(long_ref_sa_plus_consalifold_ppvs, long_ref_sa_plus_consalifold_senss, label = "Reference + ConsAlifold (Long)", marker = "^", linestyle = ":")
  line_13, = pyplot.plot(long_mafft_ginsi_plus_centroidalifold_ppvs, long_mafft_ginsi_plus_centroidalifold_senss, label = "MAFFT G-INS-i + CentroidAlifold (Long)", marker = "p", linestyle = "-")
  line_14, = pyplot.plot(long_mafft_xinsi_plus_centroidalifold_ppvs, long_mafft_xinsi_plus_centroidalifold_senss, label = "MAFFT X-INS-i + CentroidAlifold (Long)", marker = "p", linestyle = "--")
  line_15, = pyplot.plot(long_ref_sa_plus_centroidalifold_ppvs, long_ref_sa_plus_centroidalifold_senss, label = "Reference + CentroidAlifold (Long)", marker = "p", linestyle = ":")
  line_16, = pyplot.plot(long_mafft_ginsi_plus_rnaalifold_ppv, long_mafft_ginsi_plus_rnaalifold_sens, label = "MAFFT G-INS-i + RNAalifold (Long)", marker = "d", linestyle = "-")
  line_17, = pyplot.plot(long_mafft_xinsi_plus_rnaalifold_ppv, long_mafft_xinsi_plus_rnaalifold_sens, label = "MAFFT X-INS-i + RNAalifold (Long)", marker = "d", linestyle = "--")
  line_18, = pyplot.plot(long_ref_sa_plus_rnaalifold_ppv, long_ref_sa_plus_rnaalifold_sens, label = "Reference + RNAalifold (Long)", marker = "d", linestyle = ":")
  pyplot.xlabel("Positive predictive value")
  pyplot.ylabel("Sensitivity")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4, line_5, line_6, line_7, line_8, line_9], loc = 3)
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/ppvs_vs_senss_on_css_estimation.eps", bbox_inches = "tight")
  pyplot.figure()
  line_1, = pyplot.plot(short_mafft_ginsi_plus_consalifold_fprs, short_mafft_ginsi_plus_consalifold_senss, label = "MAFFT G-INS-i + ConsAlifold (Short)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(short_mafft_xinsi_plus_consalifold_fprs, short_mafft_xinsi_plus_consalifold_senss, label = "MAFFT X-INS-i + ConsAlifold (Short)", marker = "o", linestyle = "--")
  line_3, = pyplot.plot(short_ref_sa_plus_consalifold_fprs, short_ref_sa_plus_consalifold_senss, label = "Reference + ConsAlifold (Short)", marker = "o", linestyle = ":")
  line_4, = pyplot.plot(short_mafft_ginsi_plus_centroidalifold_fprs, short_mafft_ginsi_plus_centroidalifold_senss, label = "MAFFT G-INS-i + CentroidAlifold (Short)", marker = "s", linestyle = "-")
  line_5, = pyplot.plot(short_mafft_xinsi_plus_centroidalifold_fprs, short_mafft_xinsi_plus_centroidalifold_senss, label = "MAFFT X-INS-i + CentroidAlifold (Short)", marker = "s", linestyle = "--")
  line_6, = pyplot.plot(short_ref_sa_plus_centroidalifold_fprs, short_ref_sa_plus_centroidalifold_senss, label = "Reference + CentroidAlifold (Short)", marker = "s", linestyle = ":")
  line_7, = pyplot.plot(short_mafft_ginsi_plus_rnaalifold_fpr, short_mafft_ginsi_plus_rnaalifold_sens, label = "MAFFT G-INS-i + RNAAlifold (Short)", marker = "v", linestyle = "-")
  line_8, = pyplot.plot(short_mafft_xinsi_plus_rnaalifold_fpr, short_mafft_xinsi_plus_rnaalifold_sens, label = "MAFFT X-INS-i + RNAAlifold (Short)", marker = "v", linestyle = "--")
  line_9, = pyplot.plot(short_ref_sa_plus_rnaalifold_fpr, short_ref_sa_plus_rnaalifold_sens, label = "Reference + RNAAlifold (Short)", marker = "v", linestyle = ":")
  line_10, = pyplot.plot(long_mafft_ginsi_plus_consalifold_fprs, long_mafft_ginsi_plus_consalifold_senss, label = "MAFFT G-INS-i + ConsAlifold (Long)", marker = "^", linestyle = "-")
  line_11, = pyplot.plot(long_mafft_xinsi_plus_consalifold_fprs, long_mafft_xinsi_plus_consalifold_senss, label = "MAFFT X-INS-i + ConsAlifold (Long)", marker = "^", linestyle = "--")
  line_12, = pyplot.plot(long_ref_sa_plus_consalifold_fprs, long_ref_sa_plus_consalifold_senss, label = "Reference + ConsAlifold (Long)", marker = "^", linestyle = ":")
  line_13, = pyplot.plot(long_mafft_ginsi_plus_centroidalifold_fprs, long_mafft_ginsi_plus_centroidalifold_senss, label = "MAFFT G-INS-i + CentroidAlifold (Long)", marker = "p", linestyle = "-")
  line_14, = pyplot.plot(long_mafft_xinsi_plus_centroidalifold_fprs, long_mafft_xinsi_plus_centroidalifold_senss, label = "MAFFT X-INS-i + CentroidAlifold (Long)", marker = "p", linestyle = "--")
  line_15, = pyplot.plot(long_ref_sa_plus_centroidalifold_fprs, long_ref_sa_plus_centroidalifold_senss, label = "Reference + CentroidAlifold (Long)", marker = "p", linestyle = ":")
  line_16, = pyplot.plot(long_mafft_ginsi_plus_rnaalifold_fpr, long_mafft_ginsi_plus_rnaalifold_sens, label = "MAFFT G-INS-i + RNAAlifold (Long)", marker = "d", linestyle = "-")
  line_17, = pyplot.plot(long_mafft_xinsi_plus_rnaalifold_fpr, long_mafft_xinsi_plus_rnaalifold_sens, label = "MAFFT X-INS-i + RNAAlifold (Long)", marker = "d", linestyle = "--")
  line_18, = pyplot.plot(long_ref_sa_plus_rnaalifold_fpr, long_ref_sa_plus_rnaalifold_sens, label = "Reference + RNAAlifold (Long)", marker = "d", linestyle = ":")
  pyplot.xlabel("False positive rate")
  pyplot.ylabel("Sensitivity")
  pyplot.legend(handles = [line_10, line_11, line_12, line_13, line_14, line_15, line_16, line_17, line_18], loc = 4)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/fprs_vs_senss_on_css_estimation.eps", bbox_inches = "tight")

def get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css):
  num_of_rnas = len(rna_seq_lens)
  tp = fp = tn = fn = 0
  estimated_css, estimated_flat_css = estimated_css_and_flat_css
  ref_css, ref_flat_css = ref_css_and_flat_css
  for m in range(0, num_of_rnas):
    sub_estimated_css = estimated_css[m]
    sub_ref_css = ref_css[m]
    sub_estimated_flat_css = estimated_flat_css[m]
    sub_ref_flat_css = ref_flat_css[m]
    rna_seq_len_1 = rna_seq_lens[m]
    for i in range(0, rna_seq_len_1):
      estimated_bin = i in sub_estimated_flat_css
      ref_bin = i in sub_ref_flat_css
      if estimated_bin == ref_bin:
        if estimated_bin == False:
          tn += 1
      else:
        if estimated_bin == True:
          fp += 1
        else:
          fn += 1
      for j in range(i + 1, rna_seq_len_1):
        estimated_bin = (i, j) in sub_estimated_css
        ref_bin = (i, j) in sub_ref_css
        if estimated_bin == ref_bin:
          if estimated_bin == True:
            tp += 1
  return tp, tn, fp, fn

if __name__ == "__main__":
  main()
