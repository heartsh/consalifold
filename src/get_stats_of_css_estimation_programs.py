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
from math import sqrt
import multiprocessing

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  num_of_threads = multiprocessing.cpu_count()
  mafft_plus_consalifold_ppvs = []
  mafft_plus_consalifold_senss = []
  mafft_plus_consalifold_fprs = []
  probcons_plus_consalifold_ppvs = []
  probcons_plus_consalifold_senss = []
  probcons_plus_consalifold_fprs = []
  clustalw_plus_consalifold_ppvs = []
  clustalw_plus_consalifold_senss = []
  clustalw_plus_consalifold_fprs = []
  mafft_plus_centroidalifold_ppvs = []
  mafft_plus_centroidalifold_senss = []
  mafft_plus_centroidalifold_fprs = []
  probcons_plus_centroidalifold_ppvs = []
  probcons_plus_centroidalifold_senss = []
  probcons_plus_centroidalifold_fprs = []
  clustalw_plus_centroidalifold_ppvs = []
  clustalw_plus_centroidalifold_senss = []
  clustalw_plus_centroidalifold_fprs = []
  mafft_plus_rnaalifold_ppv = mafft_plus_rnaalifold_sens = mafft_plus_rnaalifold_fpr = 0.
  probcons_plus_rnaalifold_ppv = probcons_plus_rnaalifold_sens = probcons_plus_rnaalifold_fpr = 0.
  clustalw_plus_rnaalifold_ppv = clustalw_plus_rnaalifold_sens = clustalw_plus_rnaalifold_fpr = 0.
  mafft_plus_consalifold_f1_score = probcons_plus_consalifold_f1_score = clustalw_plus_consalifold_f1_score = 0.
  mafft_plus_consalifold_mcc = probcons_plus_consalifold_mcc = clustalw_plus_consalifold_mcc = 0.
  mafft_plus_centroidalifold_f1_score = probcons_plus_centroidalifold_f1_score = clustalw_plus_centroidalifold_f1_score = 0.
  mafft_plus_centroidalifold_mcc = probcons_plus_centroidalifold_mcc = clustalw_plus_centroidalifold_mcc = 0.
  mafft_plus_rnaalifold_f1_score = probcons_plus_rnaalifold_f1_score = clustalw_plus_rnaalifold_f1_score = 0.
  mafft_plus_rnaalifold_mcc = probcons_plus_rnaalifold_mcc = clustalw_plus_rnaalifold_mcc = 0.
  mafft_plus_petfold_f1_score = probcons_plus_petfold_f1_score = clustalw_plus_petfold_f1_score = 0.
  mafft_plus_petfold_mcc = probcons_plus_petfold_mcc = clustalw_plus_petfold_mcc = 0.
  gammas = [2. ** i for i in range(-7, 11)]
  rna_fam_dir_path = asset_dir_path + "/compiled_rna_fams"
  # rna_fam_dir_path = asset_dir_path + "/compiled_rna_fams_4_micro_bench"
  ref_sa_dir_path = asset_dir_path + "/ref_sas"
  # ref_sa_dir_path = asset_dir_path + "/ref_sas_4_micro_bench"
  mafft_plus_consalifold_css_dir_path = asset_dir_path + "/mafft_plus_consalifold"
  probcons_plus_consalifold_css_dir_path = asset_dir_path + "/probcons_plus_consalifold"
  clustalw_plus_consalifold_css_dir_path = asset_dir_path + "/clustalw_plus_consalifold"
  mafft_plus_centroidalifold_css_dir_path = asset_dir_path + "/mafft_plus_centroidalifold"
  probcons_plus_centroidalifold_css_dir_path = asset_dir_path + "/probcons_plus_centroidalifold"
  clustalw_plus_centroidalifold_css_dir_path = asset_dir_path + "/clustalw_plus_centroidalifold"
  mafft_plus_rnaalifold_css_dir_path = asset_dir_path + "/mafft_plus_rnaalifold"
  probcons_plus_rnaalifold_css_dir_path = asset_dir_path + "/probcons_plus_rnaalifold"
  clustalw_plus_rnaalifold_css_dir_path = asset_dir_path + "/clustalw_plus_rnaalifold"
  mafft_plus_petfold_css_dir_path = asset_dir_path + "/mafft_plus_petfold"
  probcons_plus_petfold_css_dir_path = asset_dir_path + "/probcons_plus_petfold"
  clustalw_plus_petfold_css_dir_path = asset_dir_path + "/clustalw_plus_petfold"
  pool = multiprocessing.Pool(num_of_threads)
  for gamma in gammas:
    mafft_plus_consalifold_count_params = []
    probcons_plus_consalifold_count_params = []
    clustalw_plus_consalifold_count_params = []
    mafft_plus_centroidalifold_count_params = []
    probcons_plus_centroidalifold_count_params = []
    clustalw_plus_centroidalifold_count_params = []
    mafft_plus_rnaalifold_count_params = []
    probcons_plus_rnaalifold_count_params = []
    clustalw_plus_rnaalifold_count_params = []
    mafft_plus_petfold_count_params = []
    probcons_plus_petfold_count_params = []
    clustalw_plus_petfold_count_params = []
    gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
    # mafft_plus_consalifold_tp = mafft_plus_consalifold_tn = mafft_plus_consalifold_fp = mafft_plus_consalifold_fn = 0.
    # probcons_plus_consalifold_tp = probcons_plus_consalifold_tn = probcons_plus_consalifold_fp = probcons_plus_consalifold_fn = 0.
    # clustalw_plus_consalifold_tp = clustalw_plus_consalifold_tn = clustalw_plus_consalifold_fp = clustalw_plus_consalifold_fn = 0.
    # mafft_plus_centroidalifold_tp = mafft_plus_centroidalifold_tn = mafft_plus_centroidalifold_fp = mafft_plus_centroidalifold_fn = 0.
    # probcons_plus_centroidalifold_tp = probcons_plus_centroidalifold_tn = probcons_plus_centroidalifold_fp = probcons_plus_centroidalifold_fn = 0.
    # clustalw_plus_centroidalifold_tp = clustalw_plus_centroidalifold_tn = clustalw_plus_centroidalifold_fp = clustalw_plus_centroidalifold_fn = 0.
    # mafft_plus_rnaalifold_tp = mafft_plus_rnaalifold_tn = mafft_plus_rnaalifold_fp = mafft_plus_rnaalifold_fn = 0.
    # probcons_plus_rnaalifold_tp = probcons_plus_rnaalifold_tn = probcons_plus_rnaalifold_fp = probcons_plus_rnaalifold_fn = 0.
    # clustalw_plus_rnaalifold_tp = clustalw_plus_rnaalifold_tn = clustalw_plus_rnaalifold_fp = clustalw_plus_rnaalifold_fn = 0.
    # mafft_plus_petfold_tp = mafft_plus_petfold_tn = mafft_plus_petfold_fp = mafft_plus_petfold_fn = 0.
    # probcons_plus_petfold_tp = probcons_plus_petfold_tn = probcons_plus_petfold_fp = probcons_plus_petfold_fn = 0.
    # clustalw_plus_petfold_tp = clustalw_plus_petfold_tn = clustalw_plus_petfold_fp = clustalw_plus_petfold_fn = 0.
    for rna_fam_file in os.listdir(rna_fam_dir_path):
      if not rna_fam_file.endswith(".fa"):
        continue
      rna_seq_file_path = os.path.join(rna_fam_dir_path, rna_fam_file)
      rna_seq_lens = [len(rna_seq.seq) for rna_seq in SeqIO.parse(rna_seq_file_path, "fasta")]
      num_of_rnas = len(rna_seq_lens)
      (rna_fam_name, extension) = os.path.splitext(rna_fam_file)
      ref_css_file_path = os.path.join(ref_sa_dir_path, rna_fam_name + ".sth")
      ref_css_and_flat_css = utils.get_css_and_flat_css(ref_css_file_path)
      mafft_plus_consalifold_estimated_css_dir_path = os.path.join(mafft_plus_consalifold_css_dir_path, rna_fam_name)
      probcons_plus_consalifold_estimated_css_dir_path = os.path.join(probcons_plus_consalifold_css_dir_path, rna_fam_name)
      clustalw_plus_consalifold_estimated_css_dir_path = os.path.join(clustalw_plus_consalifold_css_dir_path, rna_fam_name)
      mafft_plus_centroidalifold_estimated_css_dir_path = os.path.join(mafft_plus_centroidalifold_css_dir_path, rna_fam_name)
      probcons_plus_centroidalifold_estimated_css_dir_path = os.path.join(probcons_plus_centroidalifold_css_dir_path, rna_fam_name)
      clustalw_plus_centroidalifold_estimated_css_dir_path = os.path.join(clustalw_plus_centroidalifold_css_dir_path, rna_fam_name)
      mafft_plus_consalifold_estimated_css_file_path = os.path.join(mafft_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css_and_flat_css = utils.get_css_and_flat_css(mafft_plus_consalifold_estimated_css_file_path)
      mafft_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css))
      if False:
        tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
        mafft_plus_consalifold_tp += tp
        mafft_plus_consalifold_tn += tn
        mafft_plus_consalifold_fp += fp
        mafft_plus_consalifold_fn += fn
      probcons_plus_consalifold_estimated_css_file_path = os.path.join(probcons_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css_and_flat_css = utils.get_css_and_flat_css(probcons_plus_consalifold_estimated_css_file_path)
      probcons_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css))
      if False:
        tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
        probcons_plus_consalifold_tp += tp
        probcons_plus_consalifold_tn += tn
        probcons_plus_consalifold_fp += fp
        probcons_plus_consalifold_fn += fn
      clustalw_plus_consalifold_estimated_css_file_path = os.path.join(clustalw_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css_and_flat_css = utils.get_css_and_flat_css(clustalw_plus_consalifold_estimated_css_file_path)
      clustalw_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css))
      if False:
        tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
        clustalw_plus_consalifold_tp += tp
        clustalw_plus_consalifold_tn += tn
        clustalw_plus_consalifold_fp += fp
        clustalw_plus_consalifold_fn += fn
      mafft_plus_centroidalifold_estimated_css_file_path = os.path.join(mafft_plus_centroidalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css_and_flat_css = utils.get_css_and_flat_css(mafft_plus_centroidalifold_estimated_css_file_path)
      mafft_plus_centroidalifold_count_params.insert(0, (rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css))
      if False:
        tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
        mafft_plus_centroidalifold_tp += tp
        mafft_plus_centroidalifold_tn += tn
        mafft_plus_centroidalifold_fp += fp
        mafft_plus_centroidalifold_fn += fn
      probcons_plus_centroidalifold_estimated_css_file_path = os.path.join(probcons_plus_centroidalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css_and_flat_css = utils.get_css_and_flat_css(probcons_plus_centroidalifold_estimated_css_file_path)
      probcons_plus_centroidalifold_count_params.insert(0, (rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css))
      if False:
        tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
        probcons_plus_centroidalifold_tp += tp
        probcons_plus_centroidalifold_tn += tn
        probcons_plus_centroidalifold_fp += fp
        probcons_plus_centroidalifold_fn += fn
      clustalw_plus_centroidalifold_estimated_css_file_path = os.path.join(clustalw_plus_centroidalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css_and_flat_css = utils.get_css_and_flat_css(clustalw_plus_centroidalifold_estimated_css_file_path)
      clustalw_plus_centroidalifold_count_params.insert(0, (rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css))
      if False:
        tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
        clustalw_plus_centroidalifold_tp += tp
        clustalw_plus_centroidalifold_tn += tn
        clustalw_plus_centroidalifold_fp += fp
        clustalw_plus_centroidalifold_fn += fn
      if gamma == 1.:
        mafft_plus_rnaalifold_estimated_css_file_path = os.path.join(mafft_plus_rnaalifold_css_dir_path, rna_fam_name + ".sth")
        estimated_css_and_flat_css = utils.get_css_and_flat_css(mafft_plus_rnaalifold_estimated_css_file_path)
        mafft_plus_rnaalifold_count_params.insert(0, (rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css))
        if False:
          tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
          mafft_plus_rnaalifold_tp += tp
          mafft_plus_rnaalifold_tn += tn
          mafft_plus_rnaalifold_fp += fp
          mafft_plus_rnaalifold_fn += fn
        probcons_plus_rnaalifold_estimated_css_file_path = os.path.join(probcons_plus_rnaalifold_css_dir_path, rna_fam_name + ".sth")
        estimated_css_and_flat_css = utils.get_css_and_flat_css(probcons_plus_rnaalifold_estimated_css_file_path)
        probcons_plus_rnaalifold_count_params.insert(0, (rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css))
        if False:
          tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
          probcons_plus_rnaalifold_tp += tp
          probcons_plus_rnaalifold_tn += tn
          probcons_plus_rnaalifold_fp += fp
          probcons_plus_rnaalifold_fn += fn
        clustalw_plus_rnaalifold_estimated_css_file_path = os.path.join(clustalw_plus_rnaalifold_css_dir_path, rna_fam_name + ".sth")
        estimated_css_and_flat_css = utils.get_css_and_flat_css(clustalw_plus_rnaalifold_estimated_css_file_path)
        clustalw_plus_rnaalifold_count_params.insert(0, (rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css))
        if False:
          tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
          clustalw_plus_rnaalifold_tp += tp
          clustalw_plus_rnaalifold_tn += tn
          clustalw_plus_rnaalifold_fp += fp
          clustalw_plus_rnaalifold_fn += fn
      if gamma == 1.:
        mafft_plus_petfold_estimated_css_file_path = os.path.join(mafft_plus_petfold_css_dir_path, rna_fam_name + ".sth")
        estimated_css_and_flat_css = utils.get_css_and_flat_css(mafft_plus_petfold_estimated_css_file_path)
        mafft_plus_petfold_count_params.insert(0, (rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css))
        if False:
          tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
          mafft_plus_petfold_tp += tp
          mafft_plus_petfold_tn += tn
          mafft_plus_petfold_fp += fp
          mafft_plus_petfold_fn += fn
        probcons_plus_petfold_estimated_css_file_path = os.path.join(probcons_plus_petfold_css_dir_path, rna_fam_name + ".sth")
        estimated_css_and_flat_css = utils.get_css_and_flat_css(probcons_plus_petfold_estimated_css_file_path)
        probcons_plus_petfold_count_params.insert(0, (rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css))
        if False:
          tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
          probcons_plus_petfold_tp += tp
          probcons_plus_petfold_tn += tn
          probcons_plus_petfold_fp += fp
          probcons_plus_petfold_fn += fn
        clustalw_plus_petfold_estimated_css_file_path = os.path.join(clustalw_plus_petfold_css_dir_path, rna_fam_name + ".sth")
        estimated_css_and_flat_css = utils.get_css_and_flat_css(clustalw_plus_petfold_estimated_css_file_path)
        clustalw_plus_petfold_count_params.insert(0, (rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css))
        if False:
          tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
          clustalw_plus_petfold_tp += tp
          clustalw_plus_petfold_tn += tn
          clustalw_plus_petfold_fp += fp
          clustalw_plus_petfold_fn += fn
    results = pool.map(get_bin_counts, mafft_plus_consalifold_count_params)
    mafft_plus_consalifold_tp, mafft_plus_consalifold_tn, mafft_plus_consalifold_fp, mafft_plus_consalifold_fn = final_sum(results)
    ppv = get_ppv(mafft_plus_consalifold_tp, mafft_plus_consalifold_fp)
    sens = get_sens(mafft_plus_consalifold_tp, mafft_plus_consalifold_fn)
    fpr = get_fpr(mafft_plus_consalifold_tn, mafft_plus_consalifold_fp)
    mafft_plus_consalifold_ppvs.insert(0, ppv)
    mafft_plus_consalifold_senss.insert(0, sens)
    mafft_plus_consalifold_fprs.insert(0, fpr)
    if gamma == 1.:
      mafft_plus_consalifold_f1_score = get_f1_score(ppv, sens)
      mafft_plus_consalifold_mcc = get_mcc(mafft_plus_consalifold_tp, mafft_plus_consalifold_tn, mafft_plus_consalifold_fp, mafft_plus_consalifold_fn)
    results = pool.map(get_bin_counts, probcons_plus_consalifold_count_params)
    probcons_plus_consalifold_tp, probcons_plus_consalifold_tn, probcons_plus_consalifold_fp, probcons_plus_consalifold_fn = final_sum(results)
    ppv = get_ppv(probcons_plus_consalifold_tp, probcons_plus_consalifold_fp)
    sens = get_sens(probcons_plus_consalifold_tp, probcons_plus_consalifold_fn)
    fpr = get_fpr(probcons_plus_consalifold_tn, probcons_plus_consalifold_fp)
    probcons_plus_consalifold_ppvs.insert(0, ppv)
    probcons_plus_consalifold_senss.insert(0, sens)
    probcons_plus_consalifold_fprs.insert(0, fpr)
    if gamma == 1.:
      probcons_plus_consalifold_f1_score = get_f1_score(ppv, sens)
      probcons_plus_consalifold_mcc = get_mcc(probcons_plus_consalifold_tp, probcons_plus_consalifold_tn, probcons_plus_consalifold_fp, probcons_plus_consalifold_fn)
    results = pool.map(get_bin_counts, clustalw_plus_consalifold_count_params)
    clustalw_plus_consalifold_tp, clustalw_plus_consalifold_tn, clustalw_plus_consalifold_fp, clustalw_plus_consalifold_fn = final_sum(results)
    ppv = get_ppv(clustalw_plus_consalifold_tp, clustalw_plus_consalifold_fp)
    sens = get_sens(clustalw_plus_consalifold_tp, clustalw_plus_consalifold_fn)
    fpr = get_fpr(clustalw_plus_consalifold_tn, clustalw_plus_consalifold_fp)
    clustalw_plus_consalifold_ppvs.insert(0, ppv)
    clustalw_plus_consalifold_senss.insert(0, sens)
    clustalw_plus_consalifold_fprs.insert(0, fpr)
    if gamma == 1.:
      clustalw_plus_consalifold_f1_score = get_f1_score(ppv, sens)
      clustalw_plus_consalifold_mcc = get_mcc(clustalw_plus_consalifold_tp, clustalw_plus_consalifold_tn, clustalw_plus_consalifold_fp, clustalw_plus_consalifold_fn)
    results = pool.map(get_bin_counts, mafft_plus_centroidalifold_count_params)
    mafft_plus_centroidalifold_tp, mafft_plus_centroidalifold_tn, mafft_plus_centroidalifold_fp, mafft_plus_centroidalifold_fn = final_sum(results)
    ppv = get_ppv(mafft_plus_centroidalifold_tp, mafft_plus_centroidalifold_fp)
    sens = get_sens(mafft_plus_centroidalifold_tp, mafft_plus_centroidalifold_fn)
    fpr = get_fpr(mafft_plus_centroidalifold_tn, mafft_plus_centroidalifold_fp)
    mafft_plus_centroidalifold_ppvs.insert(0, ppv)
    mafft_plus_centroidalifold_senss.insert(0, sens)
    mafft_plus_centroidalifold_fprs.insert(0, fpr)
    if gamma == 1.:
      mafft_plus_centroidalifold_f1_score = get_f1_score(ppv, sens)
      mafft_plus_centroidalifold_mcc = get_mcc(mafft_plus_centroidalifold_tp, mafft_plus_centroidalifold_tn, mafft_plus_centroidalifold_fp, mafft_plus_centroidalifold_fn)
    results = pool.map(get_bin_counts, probcons_plus_centroidalifold_count_params)
    probcons_plus_centroidalifold_tp, probcons_plus_centroidalifold_tn, probcons_plus_centroidalifold_fp, probcons_plus_centroidalifold_fn = final_sum(results)
    ppv = get_ppv(probcons_plus_centroidalifold_tp, probcons_plus_centroidalifold_fp)
    sens = get_sens(probcons_plus_centroidalifold_tp, probcons_plus_centroidalifold_fn)
    fpr = get_fpr(probcons_plus_centroidalifold_tn, probcons_plus_centroidalifold_fp)
    probcons_plus_centroidalifold_ppvs.insert(0, ppv)
    probcons_plus_centroidalifold_senss.insert(0, sens)
    probcons_plus_centroidalifold_fprs.insert(0, fpr)
    if gamma == 1.:
      probcons_plus_centroidalifold_f1_score = get_f1_score(ppv, sens)
      probcons_plus_centroidalifold_mcc = get_mcc(probcons_plus_centroidalifold_tp, probcons_plus_centroidalifold_tn, mafft_plus_centroidalifold_fp, probcons_plus_centroidalifold_fn)
    results = pool.map(get_bin_counts, clustalw_plus_centroidalifold_count_params)
    clustalw_plus_centroidalifold_tp, clustalw_plus_centroidalifold_tn, clustalw_plus_centroidalifold_fp, clustalw_plus_centroidalifold_fn = final_sum(results)
    ppv = get_ppv(clustalw_plus_centroidalifold_tp, clustalw_plus_centroidalifold_fp)
    sens = get_sens(clustalw_plus_centroidalifold_tp, clustalw_plus_centroidalifold_fn)
    fpr = get_fpr(clustalw_plus_centroidalifold_tn, clustalw_plus_centroidalifold_fp)
    clustalw_plus_centroidalifold_ppvs.insert(0, ppv)
    clustalw_plus_centroidalifold_senss.insert(0, sens)
    clustalw_plus_centroidalifold_fprs.insert(0, fpr)
    if gamma == 1.:
      clustalw_plus_centroidalifold_f1_score = get_f1_score(ppv, sens)
      clustalw_plus_centroidalifold_mcc = get_mcc(clustalw_plus_centroidalifold_tp, clustalw_plus_centroidalifold_tn, clustalw_plus_centroidalifold_fp, clustalw_plus_centroidalifold_fn)
    if gamma == 1.:
      results = pool.map(get_bin_counts, mafft_plus_rnaalifold_count_params)
      mafft_plus_rnaalifold_tp, mafft_plus_rnaalifold_tn, mafft_plus_rnaalifold_fp, mafft_plus_rnaalifold_fn = final_sum(results)
      mafft_plus_rnaalifold_ppv = get_ppv(mafft_plus_rnaalifold_tp, mafft_plus_rnaalifold_fp)
      mafft_plus_rnaalifold_sens = get_sens(mafft_plus_rnaalifold_tp, mafft_plus_rnaalifold_fn)
      mafft_plus_rnaalifold_fpr = get_fpr(mafft_plus_rnaalifold_tn, mafft_plus_rnaalifold_fp)
      mafft_plus_rnaalifold_f1_score = get_f1_score(mafft_plus_rnaalifold_ppv, mafft_plus_rnaalifold_sens)
      mafft_plus_rnaalifold_mcc = get_mcc(mafft_plus_rnaalifold_tp, mafft_plus_rnaalifold_tn, mafft_plus_rnaalifold_fp, mafft_plus_rnaalifold_fn)
      results = pool.map(get_bin_counts, probcons_plus_rnaalifold_count_params)
      probcons_plus_rnaalifold_tp, probcons_plus_rnaalifold_tn, probcons_plus_rnaalifold_fp, probcons_plus_rnaalifold_fn = final_sum(results)
      probcons_plus_rnaalifold_ppv = get_ppv(probcons_plus_rnaalifold_tp, probcons_plus_rnaalifold_fp)
      probcons_plus_rnaalifold_sens = get_sens(probcons_plus_rnaalifold_tp, probcons_plus_rnaalifold_fn)
      probcons_plus_rnaalifold_fpr = get_fpr(probcons_plus_rnaalifold_tn, probcons_plus_rnaalifold_fp)
      probcons_plus_rnaalifold_f1_score = get_f1_score(mafft_plus_rnaalifold_ppv, mafft_plus_rnaalifold_sens)
      probcons_plus_rnaalifold_mcc = get_mcc(probcons_plus_rnaalifold_tp, probcons_plus_rnaalifold_tn, probcons_plus_rnaalifold_fp, probcons_plus_rnaalifold_fn)
      results = pool.map(get_bin_counts, clustalw_plus_rnaalifold_count_params)
      clustalw_plus_rnaalifold_tp, clustalw_plus_rnaalifold_tn, clustalw_plus_rnaalifold_fp, clustalw_plus_rnaalifold_fn = final_sum(results)
      clustalw_plus_rnaalifold_ppv = get_ppv(clustalw_plus_rnaalifold_tp, clustalw_plus_rnaalifold_fp)
      clustalw_plus_rnaalifold_sens = get_sens(clustalw_plus_rnaalifold_tp, clustalw_plus_rnaalifold_fn)
      clustalw_plus_rnaalifold_fpr = get_fpr(clustalw_plus_rnaalifold_tn, clustalw_plus_rnaalifold_fp)
      clustalw_plus_rnaalifold_f1_score = get_f1_score(clustalw_plus_rnaalifold_ppv, clustalw_plus_rnaalifold_sens)
      clustalw_plus_rnaalifold_mcc = get_mcc(clustalw_plus_rnaalifold_tp, clustalw_plus_rnaalifold_tn, clustalw_plus_rnaalifold_fp, clustalw_plus_rnaalifold_fn)
    if gamma == 1.:
      results = pool.map(get_bin_counts, mafft_plus_petfold_count_params)
      mafft_plus_petfold_tp, mafft_plus_petfold_tn, mafft_plus_petfold_fp, mafft_plus_petfold_fn = final_sum(results)
      mafft_plus_petfold_ppv = get_ppv(mafft_plus_petfold_tp, mafft_plus_petfold_fp)
      mafft_plus_petfold_sens = get_sens(mafft_plus_petfold_tp, mafft_plus_petfold_fn)
      mafft_plus_petfold_fpr = get_fpr(mafft_plus_petfold_tn, mafft_plus_petfold_fp)
      mafft_plus_petfold_f1_score = get_f1_score(mafft_plus_petfold_ppv, mafft_plus_petfold_sens)
      mafft_plus_petfold_mcc = get_mcc(mafft_plus_petfold_tp, mafft_plus_petfold_tn, mafft_plus_petfold_fp, mafft_plus_petfold_fn)
      results = pool.map(get_bin_counts, probcons_plus_petfold_count_params)
      probcons_plus_petfold_tp, probcons_plus_petfold_tn, probcons_plus_petfold_fp, probcons_plus_petfold_fn = final_sum(results)
      probcons_plus_petfold_ppv = get_ppv(probcons_plus_petfold_tp, probcons_plus_petfold_fp)
      probcons_plus_petfold_sens = get_sens(probcons_plus_petfold_tp, probcons_plus_petfold_fn)
      probcons_plus_petfold_fpr = get_fpr(probcons_plus_petfold_tn, probcons_plus_petfold_fp)
      probcons_plus_petfold_f1_score = get_f1_score(mafft_plus_petfold_ppv, mafft_plus_petfold_sens)
      probcons_plus_petfold_mcc = get_mcc(probcons_plus_petfold_tp, probcons_plus_petfold_tn, probcons_plus_petfold_fp, probcons_plus_petfold_fn)
      results = pool.map(get_bin_counts, clustalw_plus_petfold_count_params)
      clustalw_plus_petfold_tp, clustalw_plus_petfold_tn, clustalw_plus_petfold_fp, clustalw_plus_petfold_fn = final_sum(results)
      clustalw_plus_petfold_ppv = get_ppv(clustalw_plus_petfold_tp, clustalw_plus_petfold_fp)
      clustalw_plus_petfold_sens = get_sens(clustalw_plus_petfold_tp, clustalw_plus_petfold_fn)
      clustalw_plus_petfold_fpr = get_fpr(clustalw_plus_petfold_tn, clustalw_plus_petfold_fp)
      clustalw_plus_petfold_f1_score = get_f1_score(clustalw_plus_petfold_ppv, clustalw_plus_petfold_sens)
      clustalw_plus_petfold_mcc = get_mcc(clustalw_plus_petfold_tp, clustalw_plus_petfold_tn, clustalw_plus_petfold_fp, clustalw_plus_petfold_fn)
  mafft_plus_consalifold_ppvs = numpy.array(mafft_plus_consalifold_ppvs)
  mafft_plus_consalifold_senss = numpy.array(mafft_plus_consalifold_senss)
  mafft_plus_consalifold_fprs = numpy.array(mafft_plus_consalifold_fprs)
  probcons_plus_consalifold_ppvs = numpy.array(probcons_plus_consalifold_ppvs)
  probcons_plus_consalifold_senss = numpy.array(probcons_plus_consalifold_senss)
  probcons_plus_consalifold_fprs = numpy.array(probcons_plus_consalifold_fprs)
  clustalw_plus_consalifold_ppvs = numpy.array(clustalw_plus_consalifold_ppvs)
  clustalw_plus_consalifold_senss = numpy.array(clustalw_plus_consalifold_senss)
  clustalw_plus_consalifold_fprs = numpy.array(clustalw_plus_consalifold_fprs)
  mafft_plus_centroidalifold_ppvs = numpy.array(mafft_plus_centroidalifold_ppvs)
  mafft_plus_centroidalifold_senss = numpy.array(mafft_plus_centroidalifold_senss)
  mafft_plus_centroidalifold_fprs = numpy.array(mafft_plus_centroidalifold_fprs)
  probcons_plus_centroidalifold_ppvs = numpy.array(probcons_plus_centroidalifold_ppvs)
  probcons_plus_centroidalifold_senss = numpy.array(probcons_plus_centroidalifold_senss)
  probcons_plus_centroidalifold_fprs = numpy.array(probcons_plus_centroidalifold_fprs)
  clustalw_plus_centroidalifold_ppvs = numpy.array(clustalw_plus_centroidalifold_ppvs)
  clustalw_plus_centroidalifold_senss = numpy.array(clustalw_plus_centroidalifold_senss)
  clustalw_plus_centroidalifold_fprs = numpy.array(clustalw_plus_centroidalifold_fprs)
  line_1, = pyplot.plot(mafft_plus_consalifold_ppvs, mafft_plus_consalifold_senss, label = "MAFFT + ConsAlifold", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(probcons_plus_consalifold_ppvs, probcons_plus_consalifold_senss, label = "ProbCons + ConsAlifold", marker = "o", linestyle = "--")
  line_3, = pyplot.plot(clustalw_plus_consalifold_ppvs, clustalw_plus_consalifold_senss, label = "ClustalW + ConsAlifold", marker = "o", linestyle = ":")
  line_4, = pyplot.plot(mafft_plus_centroidalifold_ppvs, mafft_plus_centroidalifold_senss, label = "MAFFT + CentroidAlifold", marker = "s", linestyle = "-")
  line_5, = pyplot.plot(probcons_plus_centroidalifold_ppvs, probcons_plus_centroidalifold_senss, label = "ProbCons + CentroidAlifold", marker = "s", linestyle = "--")
  line_6, = pyplot.plot(clustalw_plus_centroidalifold_ppvs, clustalw_plus_centroidalifold_senss, label = "ClustalW + CentroidAlifold", marker = "s", linestyle = ":")
  line_7, = pyplot.plot(mafft_plus_rnaalifold_ppv, mafft_plus_rnaalifold_sens, label = "MAFFT + RNAalifold", marker = "v", linestyle = "-")
  line_8, = pyplot.plot(probcons_plus_rnaalifold_ppv, probcons_plus_rnaalifold_sens, label = "ProbCons + RNAalifold", marker = "v", linestyle = "--")
  line_9, = pyplot.plot(clustalw_plus_rnaalifold_ppv, clustalw_plus_rnaalifold_sens, label = "ClustalW + RNAalifold", marker = "v", linestyle = ":")
  line_10, = pyplot.plot(mafft_plus_petfold_ppv, mafft_plus_petfold_sens, label = "MAFFT + PETfold", marker = "^", linestyle = "-")
  line_11, = pyplot.plot(probcons_plus_petfold_ppv, probcons_plus_petfold_sens, label = "ProbCons + PETfold", marker = "^", linestyle = "--")
  line_12, = pyplot.plot(clustalw_plus_petfold_ppv, clustalw_plus_petfold_sens, label = "ClustalW + PETfold", marker = "^", linestyle = ":")
  pyplot.xlabel("Positive predictive value")
  pyplot.ylabel("Sensitivity")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4, line_5, line_6], loc = "lower left")
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/ppvs_vs_senss_on_css_estimation.eps", bbox_inches = "tight")
  pyplot.figure()
  line_1, = pyplot.plot(mafft_plus_consalifold_fprs, mafft_plus_consalifold_senss, label = "MAFFT + ConsAlifold", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(probcons_plus_consalifold_fprs, probcons_plus_consalifold_senss, label = "ProbCons + ConsAlifold", marker = "o", linestyle = "--")
  line_3, = pyplot.plot(clustalw_plus_consalifold_fprs, clustalw_plus_consalifold_senss, label = "ClustalW + ConsAlifold", marker = "o", linestyle = ":")
  line_4, = pyplot.plot(mafft_plus_centroidalifold_fprs, mafft_plus_centroidalifold_senss, label = "MAFFT + CentroidAlifold", marker = "s", linestyle = "-")
  line_5, = pyplot.plot(probcons_plus_centroidalifold_fprs, probcons_plus_centroidalifold_senss, label = "ProbCons + CentroidAlifold", marker = "s", linestyle = "--")
  line_6, = pyplot.plot(clustalw_plus_centroidalifold_fprs, clustalw_plus_centroidalifold_senss, label = "ClustalW + CentroidAlifold", marker = "s", linestyle = ":")
  line_7, = pyplot.plot(mafft_plus_rnaalifold_fpr, mafft_plus_rnaalifold_sens, label = "MAFFT + RNAAlifold", marker = "v", linestyle = "-")
  line_8, = pyplot.plot(probcons_plus_rnaalifold_fpr, probcons_plus_rnaalifold_sens, label = "ProbCons + RNAAlifold", marker = "v", linestyle = "--")
  line_9, = pyplot.plot(clustalw_plus_rnaalifold_fpr, clustalw_plus_rnaalifold_sens, label = "ClustalW + RNAAlifold", marker = "v", linestyle = ":")
  line_10, = pyplot.plot(mafft_plus_petfold_fpr, mafft_plus_petfold_sens, label = "MAFFT + PETfold", marker = "^", linestyle = "-")
  line_11, = pyplot.plot(probcons_plus_petfold_fpr, probcons_plus_petfold_sens, label = "ProbCons + PETfold", marker = "^", linestyle = "--")
  line_12, = pyplot.plot(clustalw_plus_petfold_fpr, clustalw_plus_petfold_sens, label = "ClustalW + PETfold", marker = "^", linestyle = ":")
  pyplot.xlabel("False positive rate")
  pyplot.ylabel("Sensitivity")
  pyplot.legend(handles = [line_7, line_8, line_9, line_10, line_11, line_12], loc = "lower right")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/fprs_vs_senss_on_css_estimation.eps", bbox_inches = "tight")
  print("MAFFT + ConsAlifold's MCC & F1 score = %.3f & %.3f" %(mafft_plus_consalifold_mcc, mafft_plus_consalifold_f1_score))
  print("ProbCons + ConsAlifold's MCC & F1 score = %.3f & %.3f" %(probcons_plus_consalifold_mcc, probcons_plus_consalifold_f1_score))
  print("ClustalW + ConsAlifold's MCC & F1 score = %.3f & %.3f" %(clustalw_plus_consalifold_mcc, clustalw_plus_consalifold_f1_score))
  print("MAFFT + CentroidAlifold's MCC & F1 score = %.3f & %.3f" %(mafft_plus_centroidalifold_mcc, mafft_plus_centroidalifold_f1_score))
  print("ProbCons + CentroidAlifold's MCC & F1 score = %.3f & %.3f" %(probcons_plus_centroidalifold_mcc, probcons_plus_centroidalifold_f1_score))
  print("ClustalW + CentroidAlifold's MCC & F1 score = %.3f & %.3f" %(clustalw_plus_centroidalifold_mcc, clustalw_plus_centroidalifold_f1_score))
  print("MAFFT + RNAalifold's MCC & F1 score = %.3f & %.3f" %(mafft_plus_rnaalifold_mcc, mafft_plus_rnaalifold_f1_score))
  print("ProbCons + RNAalifold's MCC & F1 score = %.3f & %.3f" %(probcons_plus_rnaalifold_mcc, probcons_plus_rnaalifold_f1_score))
  print("ClustalW + RNAalifold's MCC & F1 score = %.3f & %.3f" %(clustalw_plus_rnaalifold_mcc, clustalw_plus_rnaalifold_f1_score))
  print("MAFFT + PETfold's MCC & F1 score = %.3f & %.3f" %(mafft_plus_petfold_mcc, mafft_plus_petfold_f1_score))
  print("ProbCons + PETfold's MCC & F1 score = %.3f & %.3f" %(probcons_plus_petfold_mcc, probcons_plus_petfold_f1_score))
  print("ClustalW + PETfold's MCC & F1 score = %.3f & %.3f" %(clustalw_plus_petfold_mcc, clustalw_plus_petfold_f1_score))

def get_bin_counts(params):
  (rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css) = params
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

def final_sum(results):
  final_tp = final_tn = final_fp = final_fn = 0.
  for tp, tn, fp, fn in results:
    final_tp += tp
    final_tn += tn
    final_fp += fp
    final_fn += fn
  return (final_tp, final_tn, final_fp, final_fn)

def get_f1_score(ppv, sens):
  return 2 * ppv * sens / (ppv + sens)

def get_mcc(tp, tn, fp, fn):
  return (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

def get_ppv(tp, fp):
  return tp / (tp + fp)

def get_sens(tp, fn):
  return tp / (tp + fn)

def get_fpr(tn, fp):
  return fp / (tn + fp)

if __name__ == "__main__":
  main()
