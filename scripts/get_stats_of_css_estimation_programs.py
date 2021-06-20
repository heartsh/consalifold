#! /usr/bin/env python

import utils
from Bio import SeqIO
import seaborn
from matplotlib import pyplot
import os
import math
from math import sqrt
import multiprocessing
import numpy

seaborn.set()
pyplot.rcParams['legend.handlelength'] = 0
pyplot.rcParams['legend.fontsize'] = "large"
color_palette = seaborn.color_palette()
min_gamma = -4
max_gamma = 10
white = "#F2F2F2"
bracket_pairs = [("(", ")"), ("<", ">"), ("{", "}"), ("[", "]"), ("A", "a"), ("B", "b"), ("C", "c"), ("D", "d"), ("E", "e"), ]

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  num_of_threads = multiprocessing.cpu_count()
  mafft_plus_consalifold_ppvs = []
  mafft_plus_consalifold_senss = []
  mafft_plus_consalifold_fprs = []
  mafft_plus_consalifold_f1_scores = []
  mafft_plus_consalifold_mccs = []
  probcons_plus_consalifold_ppvs = []
  probcons_plus_consalifold_senss = []
  probcons_plus_consalifold_fprs = []
  probcons_plus_consalifold_f1_scores = []
  probcons_plus_consalifold_mccs = []
  clustalw_plus_consalifold_ppvs = []
  clustalw_plus_consalifold_senss = []
  clustalw_plus_consalifold_fprs = []
  clustalw_plus_consalifold_f1_scores = []
  clustalw_plus_consalifold_mccs = []
  mafft_xinsi_plus_consalifold_ppvs = []
  mafft_xinsi_plus_consalifold_senss = []
  mafft_xinsi_plus_consalifold_fprs = []
  mafft_xinsi_plus_consalifold_f1_scores = []
  mafft_xinsi_plus_consalifold_mccs = []
  ref_sa_plus_consalifold_ppvs = []
  ref_sa_plus_consalifold_senss = []
  ref_sa_plus_consalifold_fprs = []
  ref_sa_plus_consalifold_f1_scores = []
  ref_sa_plus_consalifold_mccs = []
  mafft_plus_centroidalifold_ppvs = []
  mafft_plus_centroidalifold_senss = []
  mafft_plus_centroidalifold_fprs = []
  mafft_plus_centroidalifold_f1_scores = []
  mafft_plus_centroidalifold_mccs = []
  probcons_plus_centroidalifold_ppvs = []
  probcons_plus_centroidalifold_senss = []
  probcons_plus_centroidalifold_fprs = []
  probcons_plus_centroidalifold_f1_scores = []
  probcons_plus_centroidalifold_mccs = []
  clustalw_plus_centroidalifold_ppvs = []
  clustalw_plus_centroidalifold_senss = []
  clustalw_plus_centroidalifold_fprs = []
  clustalw_plus_centroidalifold_f1_scores = []
  clustalw_plus_centroidalifold_mccs = []
  mafft_xinsi_plus_centroidalifold_ppvs = []
  mafft_xinsi_plus_centroidalifold_senss = []
  mafft_xinsi_plus_centroidalifold_fprs = []
  mafft_xinsi_plus_centroidalifold_f1_scores = []
  mafft_xinsi_plus_centroidalifold_mccs = []
  ref_sa_plus_centroidalifold_ppvs = []
  ref_sa_plus_centroidalifold_senss = []
  ref_sa_plus_centroidalifold_fprs = []
  ref_sa_plus_centroidalifold_f1_scores = []
  ref_sa_plus_centroidalifold_mccs = []
  mafft_plus_petfold_ppvs = []
  mafft_plus_petfold_senss = []
  mafft_plus_petfold_fprs = []
  mafft_plus_petfold_f1_scores = []
  mafft_plus_petfold_mccs = []
  probcons_plus_petfold_ppvs = []
  probcons_plus_petfold_senss = []
  probcons_plus_petfold_fprs = []
  probcons_plus_petfold_f1_scores = []
  probcons_plus_petfold_mccs = []
  clustalw_plus_petfold_ppvs = []
  clustalw_plus_petfold_senss = []
  clustalw_plus_petfold_fprs = []
  clustalw_plus_petfold_f1_scores = []
  clustalw_plus_petfold_mccs = []
  mafft_xinsi_plus_petfold_ppvs = []
  mafft_xinsi_plus_petfold_senss = []
  mafft_xinsi_plus_petfold_fprs = []
  mafft_xinsi_plus_petfold_f1_scores = []
  mafft_xinsi_plus_petfold_mccs = []
  ref_sa_plus_petfold_ppvs = []
  ref_sa_plus_petfold_senss = []
  ref_sa_plus_petfold_fprs = []
  ref_sa_plus_petfold_f1_scores = []
  ref_sa_plus_petfold_mccs = []
  mafft_plus_rnaalifold_ppv = mafft_plus_rnaalifold_sens = mafft_plus_rnaalifold_fpr = mafft_plus_rnaalifold_f1_score = mafft_plus_rnaalifold_mcc = 0.
  probcons_plus_rnaalifold_ppv = probcons_plus_rnaalifold_sens = probcons_plus_rnaalifold_fpr = probcons_plus_rnaalifold_f1_score = probcons_plus_rnaalifold_mcc = 0.
  clustalw_plus_rnaalifold_ppv = clustalw_plus_rnaalifold_sens = clustalw_plus_rnaalifold_fpr = clustalw_plus_rnaalifold_f1_score = clustalw_plus_rnaalifold_mcc = 0.
  mafft_xinsi_plus_rnaalifold_ppv = mafft_xinsi_plus_rnaalifold_sens = mafft_xinsi_plus_rnaalifold_fpr = mafft_xinsi_plus_rnaalifold_f1_score = mafft_xinsi_plus_rnaalifold_mcc = 0.
  ref_sa_plus_rnaalifold_ppv = ref_sa_plus_rnaalifold_sens = ref_sa_plus_rnaalifold_fpr = ref_sa_plus_rnaalifold_f1_score = ref_sa_plus_rnaalifold_mcc = 0.
  centroidhomfold_ppvs = []
  centroidhomfold_senss = []
  centroidhomfold_fprs = []
  centroidhomfold_f1_scores = []
  centroidhomfold_mccs = []
  locarna_ppv = locarna_sens = locarna_fpr = locarna_f1_score = locarna_mcc = 0.
  raf_ppv = raf_sens = raf_fpr = raf_f1_score = raf_mcc = 0.
  turbofold_ppvs = []
  turbofold_senss = []
  turbofold_fprs = []
  turbofold_f1_scores = []
  turbofold_mccs = []
  gammas = [2. ** i for i in range(min_gamma, max_gamma + 1)]
  rna_fam_dir_path = asset_dir_path + "/compiled_rna_fams_test"
  ref_sa_dir_path = asset_dir_path + "/ref_sas_test"
  mafft_plus_consalifold_css_dir_path = asset_dir_path + "/mafft_plus_consalifold"
  probcons_plus_consalifold_css_dir_path = asset_dir_path + "/probcons_plus_consalifold"
  clustalw_plus_consalifold_css_dir_path = asset_dir_path + "/clustalw_plus_consalifold"
  mafft_xinsi_plus_consalifold_css_dir_path = asset_dir_path + "/mafft_xinsi_plus_consalifold"
  ref_sa_plus_consalifold_css_dir_path = asset_dir_path + "/ref_sa_plus_consalifold"
  mafft_plus_centroidalifold_css_dir_path = asset_dir_path + "/mafft_plus_centroidalifold"
  probcons_plus_centroidalifold_css_dir_path = asset_dir_path + "/probcons_plus_centroidalifold"
  clustalw_plus_centroidalifold_css_dir_path = asset_dir_path + "/clustalw_plus_centroidalifold"
  mafft_xinsi_plus_centroidalifold_css_dir_path = asset_dir_path + "/mafft_xinsi_plus_centroidalifold"
  ref_sa_plus_centroidalifold_css_dir_path = asset_dir_path + "/ref_sa_plus_centroidalifold"
  mafft_plus_rnaalifold_css_dir_path = asset_dir_path + "/mafft_plus_rnaalifold"
  probcons_plus_rnaalifold_css_dir_path = asset_dir_path + "/probcons_plus_rnaalifold"
  clustalw_plus_rnaalifold_css_dir_path = asset_dir_path + "/clustalw_plus_rnaalifold"
  mafft_xinsi_plus_rnaalifold_css_dir_path = asset_dir_path + "/mafft_xinsi_plus_rnaalifold"
  ref_sa_plus_rnaalifold_css_dir_path = asset_dir_path + "/ref_sa_plus_rnaalifold"
  mafft_plus_petfold_css_dir_path = asset_dir_path + "/mafft_plus_petfold"
  probcons_plus_petfold_css_dir_path = asset_dir_path + "/probcons_plus_petfold"
  clustalw_plus_petfold_css_dir_path = asset_dir_path + "/clustalw_plus_petfold"
  mafft_xinsi_plus_petfold_css_dir_path = asset_dir_path + "/mafft_xinsi_plus_petfold"
  ref_sa_plus_petfold_css_dir_path = asset_dir_path + "/ref_sa_plus_petfold"
  centroidhomfold_ss_dir_path = asset_dir_path + "/centroidhomfold"
  locarna_css_dir_path = asset_dir_path + "/locarna"
  raf_css_dir_path = asset_dir_path + "/raf"
  turbofold_ss_dir_path = asset_dir_path + "/turbofold"
  pool = multiprocessing.Pool(num_of_threads)
  for gamma in gammas:
    mafft_plus_consalifold_count_params = []
    clustalw_plus_consalifold_count_params = []
    mafft_xinsi_plus_consalifold_count_params = []
    ref_sa_plus_consalifold_count_params = []
    probcons_plus_consalifold_count_params = []
    mafft_plus_centroidalifold_count_params = []
    probcons_plus_centroidalifold_count_params = []
    clustalw_plus_centroidalifold_count_params = []
    mafft_xinsi_plus_centroidalifold_count_params = []
    ref_sa_plus_centroidalifold_count_params = []
    mafft_plus_petfold_count_params = []
    probcons_plus_petfold_count_params = []
    clustalw_plus_petfold_count_params = []
    mafft_xinsi_plus_petfold_count_params = []
    ref_sa_plus_petfold_count_params = []
    mafft_plus_rnaalifold_count_params = []
    probcons_plus_rnaalifold_count_params = []
    clustalw_plus_rnaalifold_count_params = []
    mafft_xinsi_plus_rnaalifold_count_params = []
    ref_sa_plus_rnaalifold_count_params = []
    centroidhomfold_count_params = []
    locarna_count_params = []
    raf_count_params = []
    turbofold_count_params = []
    gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
    for rna_fam_file in os.listdir(rna_fam_dir_path):
      if not rna_fam_file.endswith(".fa"):
        continue
      rna_seq_file_path = os.path.join(rna_fam_dir_path, rna_fam_file)
      rna_seq_lens = [len(rna_seq.seq) for rna_seq in SeqIO.parse(rna_seq_file_path, "fasta")]
      num_of_rnas = len(rna_seq_lens)
      (rna_fam_name, extension) = os.path.splitext(rna_fam_file)
      ref_css_file_path = os.path.join(ref_sa_dir_path, rna_fam_name + ".sth")
      ref_css = utils.get_css(ref_css_file_path)
      mafft_plus_consalifold_estimated_css_dir_path = os.path.join(mafft_plus_consalifold_css_dir_path, rna_fam_name)
      probcons_plus_consalifold_estimated_css_dir_path = os.path.join(probcons_plus_consalifold_css_dir_path, rna_fam_name)
      clustalw_plus_consalifold_estimated_css_dir_path = os.path.join(clustalw_plus_consalifold_css_dir_path, rna_fam_name)
      mafft_xinsi_plus_consalifold_estimated_css_dir_path = os.path.join(mafft_xinsi_plus_consalifold_css_dir_path, rna_fam_name)
      ref_sa_plus_consalifold_estimated_css_dir_path = os.path.join(ref_sa_plus_consalifold_css_dir_path, rna_fam_name)
      mafft_plus_centroidalifold_estimated_css_dir_path = os.path.join(mafft_plus_centroidalifold_css_dir_path, rna_fam_name)
      probcons_plus_centroidalifold_estimated_css_dir_path = os.path.join(probcons_plus_centroidalifold_css_dir_path, rna_fam_name)
      clustalw_plus_centroidalifold_estimated_css_dir_path = os.path.join(clustalw_plus_centroidalifold_css_dir_path, rna_fam_name)
      mafft_xinsi_plus_centroidalifold_estimated_css_dir_path = os.path.join(mafft_xinsi_plus_centroidalifold_css_dir_path, rna_fam_name)
      ref_sa_plus_centroidalifold_estimated_css_dir_path = os.path.join(ref_sa_plus_centroidalifold_css_dir_path, rna_fam_name)
      mafft_plus_petfold_estimated_css_dir_path = os.path.join(mafft_plus_petfold_css_dir_path, rna_fam_name)
      probcons_plus_petfold_estimated_css_dir_path = os.path.join(probcons_plus_petfold_css_dir_path, rna_fam_name)
      clustalw_plus_petfold_estimated_css_dir_path = os.path.join(clustalw_plus_petfold_css_dir_path, rna_fam_name)
      mafft_xinsi_plus_petfold_estimated_css_dir_path = os.path.join(mafft_xinsi_plus_petfold_css_dir_path, rna_fam_name)
      ref_sa_plus_petfold_estimated_css_dir_path = os.path.join(ref_sa_plus_petfold_css_dir_path, rna_fam_name)
      centroidhomfold_estimated_ss_dir_path = os.path.join(centroidhomfold_ss_dir_path, rna_fam_name)
      turbofold_estimated_ss_dir_path = os.path.join(turbofold_ss_dir_path, rna_fam_name)
      mafft_plus_consalifold_estimated_css_file_path = os.path.join(mafft_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(mafft_plus_consalifold_estimated_css_file_path)
      mafft_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      probcons_plus_consalifold_estimated_css_file_path = os.path.join(probcons_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(probcons_plus_consalifold_estimated_css_file_path)
      probcons_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      clustalw_plus_consalifold_estimated_css_file_path = os.path.join(clustalw_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(clustalw_plus_consalifold_estimated_css_file_path)
      clustalw_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      mafft_xinsi_plus_consalifold_estimated_css_file_path = os.path.join(mafft_xinsi_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(mafft_xinsi_plus_consalifold_estimated_css_file_path)
      mafft_xinsi_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      ref_sa_plus_consalifold_estimated_css_file_path = os.path.join(ref_sa_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(ref_sa_plus_consalifold_estimated_css_file_path)
      ref_sa_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      mafft_plus_centroidalifold_estimated_css_file_path = os.path.join(mafft_plus_centroidalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(mafft_plus_centroidalifold_estimated_css_file_path)
      mafft_plus_centroidalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      probcons_plus_centroidalifold_estimated_css_file_path = os.path.join(probcons_plus_centroidalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(probcons_plus_centroidalifold_estimated_css_file_path)
      probcons_plus_centroidalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      clustalw_plus_centroidalifold_estimated_css_file_path = os.path.join(clustalw_plus_centroidalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(clustalw_plus_centroidalifold_estimated_css_file_path)
      clustalw_plus_centroidalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      mafft_xinsi_plus_centroidalifold_estimated_css_file_path = os.path.join(mafft_xinsi_plus_centroidalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(mafft_xinsi_plus_centroidalifold_estimated_css_file_path)
      mafft_xinsi_plus_centroidalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      ref_sa_plus_centroidalifold_estimated_css_file_path = os.path.join(ref_sa_plus_centroidalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(ref_sa_plus_centroidalifold_estimated_css_file_path)
      ref_sa_plus_centroidalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      mafft_plus_petfold_estimated_css_file_path = os.path.join(mafft_plus_petfold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(mafft_plus_petfold_estimated_css_file_path)
      mafft_plus_petfold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      probcons_plus_petfold_estimated_css_file_path = os.path.join(probcons_plus_petfold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(probcons_plus_petfold_estimated_css_file_path)
      probcons_plus_petfold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      clustalw_plus_petfold_estimated_css_file_path = os.path.join(clustalw_plus_petfold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(clustalw_plus_petfold_estimated_css_file_path)
      clustalw_plus_petfold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      mafft_xinsi_plus_petfold_estimated_css_file_path = os.path.join(mafft_xinsi_plus_petfold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(mafft_xinsi_plus_petfold_estimated_css_file_path)
      mafft_xinsi_plus_petfold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      ref_sa_plus_petfold_estimated_css_file_path = os.path.join(ref_sa_plus_petfold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(ref_sa_plus_petfold_estimated_css_file_path)
      ref_sa_plus_petfold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      centroidhomfold_estimated_ss_file_path = os.path.join(centroidhomfold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      estimated_sss = get_sss(centroidhomfold_estimated_ss_file_path)
      centroidhomfold_count_params.insert(0, (rna_seq_lens, estimated_sss, ref_css))
      if gamma == 1.:
        mafft_plus_rnaalifold_estimated_css_file_path = os.path.join(mafft_plus_rnaalifold_css_dir_path, rna_fam_name + ".sth")
        estimated_css = utils.get_css(mafft_plus_rnaalifold_estimated_css_file_path)
        mafft_plus_rnaalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
        probcons_plus_rnaalifold_estimated_css_file_path = os.path.join(probcons_plus_rnaalifold_css_dir_path, rna_fam_name + ".sth")
        estimated_css = utils.get_css(probcons_plus_rnaalifold_estimated_css_file_path)
        probcons_plus_rnaalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
        clustalw_plus_rnaalifold_estimated_css_file_path = os.path.join(clustalw_plus_rnaalifold_css_dir_path, rna_fam_name + ".sth")
        estimated_css = utils.get_css(clustalw_plus_rnaalifold_estimated_css_file_path)
        clustalw_plus_rnaalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
        mafft_xinsi_plus_rnaalifold_estimated_css_file_path = os.path.join(mafft_xinsi_plus_rnaalifold_css_dir_path, rna_fam_name + ".sth")
        estimated_css = utils.get_css(mafft_xinsi_plus_rnaalifold_estimated_css_file_path)
        mafft_xinsi_plus_rnaalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
        ref_sa_plus_rnaalifold_estimated_css_file_path = os.path.join(ref_sa_plus_rnaalifold_css_dir_path, rna_fam_name + ".sth")
        estimated_css = utils.get_css(ref_sa_plus_rnaalifold_estimated_css_file_path)
        ref_sa_plus_rnaalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
        locarna_estimated_css_file_path = os.path.join(locarna_css_dir_path, rna_fam_name + "/results/result.stk")
        estimated_css = utils.get_css(locarna_estimated_css_file_path)
        locarna_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
        raf_estimated_css_file_path = os.path.join(raf_css_dir_path, rna_fam_name + ".sth")
        estimated_css = utils.get_css(raf_estimated_css_file_path)
        raf_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      turbofold_estimated_ss_file_path = os.path.join(turbofold_estimated_ss_dir_path, "gamma=" + gamma_str + ".fa")
      estimated_sss = get_sss(turbofold_estimated_ss_file_path)
      turbofold_count_params.insert(0, (rna_seq_lens, estimated_sss, ref_css))
    results = pool.map(get_bin_counts, mafft_plus_consalifold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    mafft_plus_consalifold_ppvs.insert(0, ppv)
    mafft_plus_consalifold_senss.insert(0, sens)
    mafft_plus_consalifold_fprs.insert(0, fpr)
    mafft_plus_consalifold_f1_scores.append(f1_score)
    mafft_plus_consalifold_mccs.append(mcc)
    results = pool.map(get_bin_counts, probcons_plus_consalifold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    probcons_plus_consalifold_ppvs.insert(0, ppv)
    probcons_plus_consalifold_senss.insert(0, sens)
    probcons_plus_consalifold_fprs.insert(0, fpr)
    probcons_plus_consalifold_f1_scores.append(f1_score)
    probcons_plus_consalifold_mccs.append(mcc)
    results = pool.map(get_bin_counts, clustalw_plus_consalifold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    clustalw_plus_consalifold_ppvs.insert(0, ppv)
    clustalw_plus_consalifold_senss.insert(0, sens)
    clustalw_plus_consalifold_fprs.insert(0, fpr)
    clustalw_plus_consalifold_f1_scores.append(f1_score)
    clustalw_plus_consalifold_mccs.append(mcc)
    results = pool.map(get_bin_counts, mafft_xinsi_plus_consalifold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    mafft_xinsi_plus_consalifold_ppvs.insert(0, ppv)
    mafft_xinsi_plus_consalifold_senss.insert(0, sens)
    mafft_xinsi_plus_consalifold_fprs.insert(0, fpr)
    mafft_xinsi_plus_consalifold_f1_scores.append(f1_score)
    mafft_xinsi_plus_consalifold_mccs.append(mcc)
    results = pool.map(get_bin_counts, ref_sa_plus_consalifold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    ref_sa_plus_consalifold_ppvs.insert(0, ppv)
    ref_sa_plus_consalifold_senss.insert(0, sens)
    ref_sa_plus_consalifold_fprs.insert(0, fpr)
    ref_sa_plus_consalifold_f1_scores.append(f1_score)
    ref_sa_plus_consalifold_mccs.append(mcc)
    results = pool.map(get_bin_counts, mafft_plus_centroidalifold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    mafft_plus_centroidalifold_ppvs.insert(0, ppv)
    mafft_plus_centroidalifold_senss.insert(0, sens)
    mafft_plus_centroidalifold_fprs.insert(0, fpr)
    mafft_plus_centroidalifold_f1_scores.append(f1_score)
    mafft_plus_centroidalifold_mccs.append(mcc)
    results = pool.map(get_bin_counts, probcons_plus_centroidalifold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    probcons_plus_centroidalifold_ppvs.insert(0, ppv)
    probcons_plus_centroidalifold_senss.insert(0, sens)
    probcons_plus_centroidalifold_fprs.insert(0, fpr)
    probcons_plus_centroidalifold_f1_scores.append(f1_score)
    probcons_plus_centroidalifold_mccs.append(mcc)
    results = pool.map(get_bin_counts, clustalw_plus_centroidalifold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    clustalw_plus_centroidalifold_ppvs.insert(0, ppv)
    clustalw_plus_centroidalifold_senss.insert(0, sens)
    clustalw_plus_centroidalifold_fprs.insert(0, fpr)
    clustalw_plus_centroidalifold_f1_scores.append(f1_score)
    clustalw_plus_centroidalifold_mccs.append(mcc)
    results = pool.map(get_bin_counts, mafft_xinsi_plus_centroidalifold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    mafft_xinsi_plus_centroidalifold_ppvs.insert(0, ppv)
    mafft_xinsi_plus_centroidalifold_senss.insert(0, sens)
    mafft_xinsi_plus_centroidalifold_fprs.insert(0, fpr)
    mafft_xinsi_plus_centroidalifold_f1_scores.append(f1_score)
    mafft_xinsi_plus_centroidalifold_mccs.append(mcc)
    results = pool.map(get_bin_counts, ref_sa_plus_centroidalifold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    ref_sa_plus_centroidalifold_ppvs.insert(0, ppv)
    ref_sa_plus_centroidalifold_senss.insert(0, sens)
    ref_sa_plus_centroidalifold_fprs.insert(0, fpr)
    ref_sa_plus_centroidalifold_f1_scores.append(f1_score)
    ref_sa_plus_centroidalifold_mccs.append(mcc)
    results = pool.map(get_bin_counts, mafft_plus_petfold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    mafft_plus_petfold_ppvs.insert(0, ppv)
    mafft_plus_petfold_senss.insert(0, sens)
    mafft_plus_petfold_fprs.insert(0, fpr)
    mafft_plus_petfold_f1_scores.append(f1_score)
    mafft_plus_petfold_mccs.append(mcc)
    results = pool.map(get_bin_counts, probcons_plus_petfold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    probcons_plus_petfold_ppvs.insert(0, ppv)
    probcons_plus_petfold_senss.insert(0, sens)
    probcons_plus_petfold_fprs.insert(0, fpr)
    probcons_plus_petfold_f1_scores.append(f1_score)
    probcons_plus_petfold_mccs.append(mcc)
    results = pool.map(get_bin_counts, clustalw_plus_petfold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    clustalw_plus_petfold_ppvs.insert(0, ppv)
    clustalw_plus_petfold_senss.insert(0, sens)
    clustalw_plus_petfold_fprs.insert(0, fpr)
    clustalw_plus_petfold_f1_scores.append(f1_score)
    clustalw_plus_petfold_mccs.append(mcc)
    results = pool.map(get_bin_counts, mafft_xinsi_plus_petfold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    mafft_xinsi_plus_petfold_ppvs.insert(0, ppv)
    mafft_xinsi_plus_petfold_senss.insert(0, sens)
    mafft_xinsi_plus_petfold_fprs.insert(0, fpr)
    mafft_xinsi_plus_petfold_f1_scores.append(f1_score)
    mafft_xinsi_plus_petfold_mccs.append(mcc)
    results = pool.map(get_bin_counts, ref_sa_plus_petfold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    ref_sa_plus_petfold_ppvs.insert(0, ppv)
    ref_sa_plus_petfold_senss.insert(0, sens)
    ref_sa_plus_petfold_fprs.insert(0, fpr)
    ref_sa_plus_petfold_f1_scores.append(f1_score)
    ref_sa_plus_petfold_mccs.append(mcc)
    results = pool.map(get_bin_counts, centroidhomfold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    centroidhomfold_ppvs.insert(0, ppv)
    centroidhomfold_senss.insert(0, sens)
    centroidhomfold_fprs.insert(0, fpr)
    centroidhomfold_f1_scores.append(f1_score)
    centroidhomfold_mccs.append(mcc)
    if gamma == 1.:
      results = pool.map(get_bin_counts, mafft_plus_rnaalifold_count_params)
      mafft_plus_rnaalifold_ppv, mafft_plus_rnaalifold_sens, mafft_plus_rnaalifold_fpr, mafft_plus_rnaalifold_f1_score, mafft_plus_rnaalifold_mcc = get_metrics(final_sum(results))
      results = pool.map(get_bin_counts, probcons_plus_rnaalifold_count_params)
      probcons_plus_rnaalifold_ppv, probcons_plus_rnaalifold_sens, probcons_plus_rnaalifold_fpr, probcons_plus_rnaalifold_f1_score, probcons_plus_rnaalifold_mcc = get_metrics(final_sum(results))
      results = pool.map(get_bin_counts, clustalw_plus_rnaalifold_count_params)
      clustalw_plus_rnaalifold_ppv, clustalw_plus_rnaalifold_sens, clustalw_plus_rnaalifold_fpr, clustalw_plus_rnaalifold_f1_score, clustalw_plus_rnaalifold_mcc = get_metrics(final_sum(results))
      results = pool.map(get_bin_counts, mafft_xinsi_plus_rnaalifold_count_params)
      mafft_xinsi_plus_rnaalifold_ppv, mafft_xinsi_plus_rnaalifold_sens, mafft_xinsi_plus_rnaalifold_fpr, mafft_xinsi_plus_rnaalifold_f1_score, mafft_xinsi_plus_rnaalifold_mcc = get_metrics(final_sum(results))
      results = pool.map(get_bin_counts, ref_sa_plus_rnaalifold_count_params)
      ref_sa_plus_rnaalifold_ppv, ref_sa_plus_rnaalifold_sens, ref_sa_plus_rnaalifold_fpr, ref_sa_plus_rnaalifold_f1_score, ref_sa_plus_rnaalifold_mcc = get_metrics(final_sum(results))
      results = pool.map(get_bin_counts, locarna_count_params)
      locarna_ppv, locarna_sens, locarna_fpr, locarna_f1_score, locarna_mcc = get_metrics(final_sum(results))
      results = pool.map(get_bin_counts, raf_count_params)
      raf_ppv, raf_sens, raf_fpr, raf_f1_score, raf_mcc = get_metrics(final_sum(results))
    results = pool.map(get_bin_counts, turbofold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    turbofold_ppvs.insert(0, ppv)
    turbofold_senss.insert(0, sens)
    turbofold_fprs.insert(0, fpr)
    turbofold_f1_scores.append(f1_score)
    turbofold_mccs.append(mcc)
  # Figure for ProbCons.
  line_1, = pyplot.plot(probcons_plus_consalifold_ppvs, probcons_plus_consalifold_senss, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(probcons_plus_centroidalifold_ppvs, probcons_plus_centroidalifold_senss, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(probcons_plus_petfold_ppvs, probcons_plus_petfold_senss, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(probcons_plus_rnaalifold_ppv, probcons_plus_rnaalifold_sens, label = "RNAalifold", marker = "v", linestyle = ":")
  pyplot.xlabel("Precision")
  pyplot.ylabel("Recall")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4], loc = "lower left")
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/pr_curves_on_css_estimation_probcons.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for MAFFT.
  pyplot.figure()
  line_1, = pyplot.plot(mafft_plus_consalifold_ppvs, mafft_plus_consalifold_senss, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(mafft_plus_centroidalifold_ppvs, mafft_plus_centroidalifold_senss, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(mafft_plus_petfold_ppvs, mafft_plus_petfold_senss, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(mafft_plus_rnaalifold_ppv, mafft_plus_rnaalifold_sens, label = "RNAalifold", marker = "v", linestyle = ":")
  pyplot.xlabel("Precision")
  pyplot.ylabel("Recall")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4], loc = "lower left")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/pr_curves_on_css_estimation_mafft.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for ClustalW.
  pyplot.figure()
  line_1, = pyplot.plot(clustalw_plus_consalifold_ppvs, clustalw_plus_consalifold_senss, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(clustalw_plus_centroidalifold_ppvs, clustalw_plus_centroidalifold_senss, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(clustalw_plus_petfold_ppvs, clustalw_plus_petfold_senss, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(clustalw_plus_rnaalifold_ppv, clustalw_plus_rnaalifold_sens, label = "RNAalifold", marker = "v", linestyle = ":")
  pyplot.xlabel("Precision")
  pyplot.ylabel("Recall")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4], loc = "lower left")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/pr_curves_on_css_estimation_clustalw.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for MAFFT X-INS-i.
  pyplot.figure()
  line_1, = pyplot.plot(mafft_xinsi_plus_consalifold_ppvs, mafft_xinsi_plus_consalifold_senss, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(mafft_xinsi_plus_centroidalifold_ppvs, mafft_xinsi_plus_centroidalifold_senss, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(mafft_xinsi_plus_petfold_ppvs, mafft_xinsi_plus_petfold_senss, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(mafft_xinsi_plus_rnaalifold_ppv, mafft_xinsi_plus_rnaalifold_sens, label = "RNAalifold", marker = "v", linestyle = ":")
  pyplot.xlabel("Precision")
  pyplot.ylabel("Recall")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4], loc = "lower left")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/pr_curves_on_css_estimation_mafft_xinsi.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for reference sequence alignments.
  pyplot.figure()
  line_1, = pyplot.plot(ref_sa_plus_consalifold_ppvs, ref_sa_plus_consalifold_senss, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(ref_sa_plus_centroidalifold_ppvs, ref_sa_plus_centroidalifold_senss, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(ref_sa_plus_petfold_ppvs, ref_sa_plus_petfold_senss, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(ref_sa_plus_rnaalifold_ppv, ref_sa_plus_rnaalifold_sens, label = "RNAalifold", marker = "v", linestyle = ":")
  pyplot.xlabel("Precision")
  pyplot.ylabel("Recall")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4], loc = "lower left")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/pr_curves_on_css_estimation_ref_sa.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for ProbCons.
  pyplot.figure()
  line_1, = pyplot.plot(probcons_plus_consalifold_fprs, probcons_plus_consalifold_senss, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(probcons_plus_centroidalifold_fprs, probcons_plus_centroidalifold_senss, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(probcons_plus_petfold_fprs, probcons_plus_petfold_senss, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(probcons_plus_rnaalifold_fpr, probcons_plus_rnaalifold_sens, label = "RNAalifold", marker = "v", linestyle = ":")
  pyplot.xlabel("Fall-out")
  pyplot.ylabel("Recall")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/roc_curves_on_css_estimation_probcons.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for MAFFT.
  pyplot.figure()
  line_1, = pyplot.plot(mafft_plus_consalifold_fprs, mafft_plus_consalifold_senss, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(mafft_plus_centroidalifold_fprs, mafft_plus_centroidalifold_senss, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(mafft_plus_petfold_fprs, mafft_plus_petfold_senss, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(mafft_plus_rnaalifold_fpr, mafft_plus_rnaalifold_sens, label = "RNAalifold", marker = "v", linestyle = ":")
  pyplot.xlabel("Fall-out")
  pyplot.ylabel("Recall")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/roc_curves_on_css_estimation_mafft.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for ClustalW.
  pyplot.figure()
  line_1, = pyplot.plot(clustalw_plus_consalifold_fprs, clustalw_plus_consalifold_senss, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(clustalw_plus_centroidalifold_fprs, clustalw_plus_centroidalifold_senss, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(clustalw_plus_petfold_fprs, clustalw_plus_petfold_senss, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(clustalw_plus_rnaalifold_fpr, clustalw_plus_rnaalifold_sens, label = "RNAalifold", marker = "v", linestyle = ":")
  pyplot.xlabel("Fall-out")
  pyplot.ylabel("Recall")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/roc_curves_on_css_estimation_clustalw.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for MAFFT X-INS-i.
  pyplot.figure()
  line_1, = pyplot.plot(mafft_xinsi_plus_consalifold_fprs, mafft_xinsi_plus_consalifold_senss, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(mafft_xinsi_plus_centroidalifold_fprs, mafft_xinsi_plus_centroidalifold_senss, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(mafft_xinsi_plus_petfold_fprs, mafft_xinsi_plus_petfold_senss, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(mafft_xinsi_plus_rnaalifold_fpr, mafft_xinsi_plus_rnaalifold_sens, label = "RNAalifold", marker = "v", linestyle = ":")
  pyplot.xlabel("Fall-out")
  pyplot.ylabel("Recall")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/roc_curves_on_css_estimation_mafft_xinsi.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for reference sequence alignments.
  pyplot.figure()
  line_1, = pyplot.plot(ref_sa_plus_consalifold_fprs, ref_sa_plus_consalifold_senss, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(ref_sa_plus_centroidalifold_fprs, ref_sa_plus_centroidalifold_senss, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(ref_sa_plus_petfold_fprs, ref_sa_plus_petfold_senss, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(ref_sa_plus_rnaalifold_fpr, ref_sa_plus_rnaalifold_sens, label = "RNAalifold", marker = "v", linestyle = ":")
  pyplot.xlabel("Fall-out")
  pyplot.ylabel("Recall")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/roc_curves_on_css_estimation_ref_sa.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for ProbCons.
  gammas = [i for i in range(min_gamma, max_gamma + 1)]
  pyplot.figure()
  line_1, = pyplot.plot(gammas, probcons_plus_consalifold_f1_scores, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, probcons_plus_centroidalifold_f1_scores, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(gammas, probcons_plus_petfold_f1_scores, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(-2., probcons_plus_rnaalifold_f1_score, label = "RNAalifold", marker = "v", linestyle = ":")
  line_5, = pyplot.plot(min_gamma + numpy.argmax(probcons_plus_consalifold_f1_scores), max(probcons_plus_consalifold_f1_scores), label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[0])
  line_6, = pyplot.plot(min_gamma + numpy.argmax(probcons_plus_centroidalifold_f1_scores), max(probcons_plus_centroidalifold_f1_scores), label = "CentroidAlifold", marker = "s", linestyle = "--", markerfacecolor = white, markeredgecolor = color_palette[1])
  line_7, = pyplot.plot(min_gamma + numpy.argmax(probcons_plus_petfold_f1_scores), max(probcons_plus_petfold_f1_scores), label = "PETfold", marker = "^", linestyle = "-.", markerfacecolor = white, markeredgecolor = color_palette[2])
  pyplot.legend(handles = [line_5, line_6, line_7], loc = "lower right")
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("F1 score")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_f1_scores_on_css_estimation_probcons.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for MAFFT.
  pyplot.figure()
  line_1, = pyplot.plot(gammas, mafft_plus_consalifold_f1_scores, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, mafft_plus_centroidalifold_f1_scores, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(gammas, mafft_plus_petfold_f1_scores, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(-2., mafft_plus_rnaalifold_f1_score, label = "RNAalifold", marker = "v", linestyle = ":")
  line_5, = pyplot.plot(min_gamma + numpy.argmax(mafft_plus_consalifold_f1_scores), max(mafft_plus_consalifold_f1_scores), label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[0])
  line_6, = pyplot.plot(min_gamma + numpy.argmax(mafft_plus_centroidalifold_f1_scores), max(mafft_plus_centroidalifold_f1_scores), label = "CentroidAlifold", marker = "s", linestyle = "--", markerfacecolor = white, markeredgecolor = color_palette[1])
  line_7, = pyplot.plot(min_gamma + numpy.argmax(mafft_plus_petfold_f1_scores), max(mafft_plus_petfold_f1_scores), label = "PETfold", marker = "^", linestyle = "-.", markerfacecolor = white, markeredgecolor = color_palette[2])
  pyplot.legend(handles = [line_5, line_6, line_7], loc = "lower right")
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("F1 score")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_f1_scores_on_css_estimation_mafft.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for ClustalW.
  pyplot.figure()
  line_1, = pyplot.plot(gammas, clustalw_plus_consalifold_f1_scores, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, clustalw_plus_centroidalifold_f1_scores, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(gammas, clustalw_plus_petfold_f1_scores, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(-2., clustalw_plus_rnaalifold_f1_score, label = "RNAalifold", marker = "v", linestyle = ":")
  line_5, = pyplot.plot(min_gamma + numpy.argmax(clustalw_plus_consalifold_f1_scores), max(clustalw_plus_consalifold_f1_scores), label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[0])
  line_6, = pyplot.plot(min_gamma + numpy.argmax(clustalw_plus_centroidalifold_f1_scores), max(clustalw_plus_centroidalifold_f1_scores), label = "CentroidAlifold", marker = "s", linestyle = "--", markerfacecolor = white, markeredgecolor = color_palette[1])
  line_7, = pyplot.plot(min_gamma + numpy.argmax(clustalw_plus_petfold_f1_scores), max(clustalw_plus_petfold_f1_scores), label = "PETfold", marker = "^", linestyle = "-.", markerfacecolor = white, markeredgecolor = color_palette[2])
  pyplot.legend(handles = [line_5, line_6, line_7], loc = "lower right")
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("F1 score")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_f1_scores_on_css_estimation_clustalw.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for MAFFT X-INS-i.
  pyplot.figure()
  line_1, = pyplot.plot(gammas, mafft_xinsi_plus_consalifold_f1_scores, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, mafft_xinsi_plus_centroidalifold_f1_scores, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(gammas, mafft_xinsi_plus_petfold_f1_scores, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(-2., mafft_xinsi_plus_rnaalifold_f1_score, label = "RNAalifold", marker = "v", linestyle = ":")
  line_5, = pyplot.plot(min_gamma + numpy.argmax(mafft_xinsi_plus_consalifold_f1_scores), max(mafft_xinsi_plus_consalifold_f1_scores), label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[0])
  line_6, = pyplot.plot(min_gamma + numpy.argmax(mafft_xinsi_plus_centroidalifold_f1_scores), max(mafft_xinsi_plus_centroidalifold_f1_scores), label = "CentroidAlifold", marker = "s", linestyle = "--", markerfacecolor = white, markeredgecolor = color_palette[1])
  line_7, = pyplot.plot(min_gamma + numpy.argmax(mafft_xinsi_plus_petfold_f1_scores), max(mafft_xinsi_plus_petfold_f1_scores), label = "PETfold", marker = "^", linestyle = "-.", markerfacecolor = white, markeredgecolor = color_palette[2])
  pyplot.legend(handles = [line_5, line_6, line_7], loc = "lower right")
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("F1 score")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_f1_scores_on_css_estimation_mafft_xinsi.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for reference sequence alignments.
  pyplot.figure()
  line_1, = pyplot.plot(gammas, ref_sa_plus_consalifold_f1_scores, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, ref_sa_plus_centroidalifold_f1_scores, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(gammas, ref_sa_plus_petfold_f1_scores, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(-2., ref_sa_plus_rnaalifold_f1_score, label = "RNAalifold", marker = "v", linestyle = ":")
  line_5, = pyplot.plot(min_gamma + numpy.argmax(ref_sa_plus_consalifold_f1_scores), max(ref_sa_plus_consalifold_f1_scores), label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[0])
  line_6, = pyplot.plot(min_gamma + numpy.argmax(ref_sa_plus_centroidalifold_f1_scores), max(ref_sa_plus_centroidalifold_f1_scores), label = "CentroidAlifold", marker = "s", linestyle = "--", markerfacecolor = white, markeredgecolor = color_palette[1])
  line_7, = pyplot.plot(min_gamma + numpy.argmax(ref_sa_plus_petfold_f1_scores), max(ref_sa_plus_petfold_f1_scores), label = "PETfold", marker = "^", linestyle = "-.", markerfacecolor = white, markeredgecolor = color_palette[2])
  pyplot.legend(handles = [line_5, line_6, line_7], loc = "lower right")
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("F1 score")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_f1_scores_on_css_estimation_ref_sa.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for ProbCons.
  pyplot.figure()
  line_1, = pyplot.plot(gammas, probcons_plus_consalifold_mccs, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, probcons_plus_centroidalifold_mccs, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(gammas, probcons_plus_petfold_mccs, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(-2., probcons_plus_rnaalifold_mcc, label = "RNAalifold", marker = "v", linestyle = ":")
  line_5, = pyplot.plot(min_gamma + numpy.argmax(probcons_plus_consalifold_mccs), max(probcons_plus_consalifold_mccs), label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[0])
  line_6, = pyplot.plot(min_gamma + numpy.argmax(probcons_plus_centroidalifold_mccs), max(probcons_plus_centroidalifold_mccs), label = "CentroidAlifold", marker = "s", linestyle = "--", markerfacecolor = white, markeredgecolor = color_palette[1])
  line_7, = pyplot.plot(min_gamma + numpy.argmax(probcons_plus_petfold_mccs), max(probcons_plus_petfold_mccs), label = "PETfold", marker = "^", linestyle = "-.", markerfacecolor = white, markeredgecolor = color_palette[2])
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("Matthews correlation coefficient")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_mccs_on_css_estimation_probcons.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for MAFFT.
  pyplot.figure()
  line_1, = pyplot.plot(gammas, mafft_plus_consalifold_mccs, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, mafft_plus_centroidalifold_mccs, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(gammas, mafft_plus_petfold_mccs, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(-2., mafft_plus_rnaalifold_mcc, label = "RNAalifold", marker = "v", linestyle = ":")
  line_5, = pyplot.plot(min_gamma + numpy.argmax(mafft_plus_consalifold_mccs), max(mafft_plus_consalifold_mccs), label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[0])
  line_6, = pyplot.plot(min_gamma + numpy.argmax(mafft_plus_centroidalifold_mccs), max(mafft_plus_centroidalifold_mccs), label = "CentroidAlifold", marker = "s", linestyle = "--", markerfacecolor = white, markeredgecolor = color_palette[1])
  line_7, = pyplot.plot(min_gamma + numpy.argmax(mafft_plus_petfold_mccs), max(mafft_plus_petfold_mccs), label = "PETfold", marker = "^", linestyle = "-.", markerfacecolor = white, markeredgecolor = color_palette[2])
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("Matthews correlation coefficient")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_mccs_on_css_estimation_mafft.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for ClustalW.
  pyplot.figure()
  line_1, = pyplot.plot(gammas, clustalw_plus_consalifold_mccs, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, clustalw_plus_centroidalifold_mccs, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(gammas, clustalw_plus_petfold_mccs, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(-2., clustalw_plus_rnaalifold_mcc, label = "RNAalifold", marker = "v", linestyle = ":")
  line_5, = pyplot.plot(min_gamma + numpy.argmax(clustalw_plus_consalifold_mccs), max(clustalw_plus_consalifold_mccs), label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[0])
  line_6, = pyplot.plot(min_gamma + numpy.argmax(clustalw_plus_centroidalifold_mccs), max(clustalw_plus_centroidalifold_mccs), label = "CentroidAlifold", marker = "s", linestyle = "--", markerfacecolor = white, markeredgecolor = color_palette[1])
  line_7, = pyplot.plot(min_gamma + numpy.argmax(clustalw_plus_petfold_mccs), max(clustalw_plus_petfold_mccs), label = "PETfold", marker = "^", linestyle = "-.", markerfacecolor = white, markeredgecolor = color_palette[2])
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("Matthews correlation coefficient")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_mccs_on_css_estimation_clustalw.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for MAFFT X-INS-i.
  pyplot.figure()
  line_1, = pyplot.plot(gammas, mafft_xinsi_plus_consalifold_mccs, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, mafft_xinsi_plus_centroidalifold_mccs, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(gammas, mafft_xinsi_plus_petfold_mccs, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(-2., mafft_xinsi_plus_rnaalifold_mcc, label = "RNAalifold", marker = "v", linestyle = ":")
  line_5, = pyplot.plot(min_gamma + numpy.argmax(mafft_xinsi_plus_consalifold_mccs), max(mafft_xinsi_plus_consalifold_mccs), label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[0])
  line_6, = pyplot.plot(min_gamma + numpy.argmax(mafft_xinsi_plus_centroidalifold_mccs), max(mafft_xinsi_plus_centroidalifold_mccs), label = "CentroidAlifold", marker = "s", linestyle = "--", markerfacecolor = white, markeredgecolor = color_palette[1])
  line_7, = pyplot.plot(min_gamma + numpy.argmax(mafft_xinsi_plus_petfold_mccs), max(mafft_xinsi_plus_petfold_mccs), label = "PETfold", marker = "^", linestyle = "-.", markerfacecolor = white, markeredgecolor = color_palette[2])
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("Matthews correlation coefficient")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_mccs_on_css_estimation_mafft_xinsi.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for reference sequence alignments.
  pyplot.figure()
  line_1, = pyplot.plot(gammas, ref_sa_plus_consalifold_mccs, label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, ref_sa_plus_centroidalifold_mccs, label = "CentroidAlifold", marker = "s", linestyle = "--")
  line_3, = pyplot.plot(gammas, ref_sa_plus_petfold_mccs, label = "PETfold", marker = "^", linestyle = "-.")
  line_4, = pyplot.plot(-2., ref_sa_plus_rnaalifold_mcc, label = "RNAalifold", marker = "v", linestyle = ":")
  line_5, = pyplot.plot(min_gamma + numpy.argmax(ref_sa_plus_consalifold_mccs), max(ref_sa_plus_consalifold_mccs), label = "ConsAlifold (ConsProb)", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[0])
  line_6, = pyplot.plot(min_gamma + numpy.argmax(ref_sa_plus_centroidalifold_mccs), max(ref_sa_plus_centroidalifold_mccs), label = "CentroidAlifold", marker = "s", linestyle = "--", markerfacecolor = white, markeredgecolor = color_palette[1])
  line_7, = pyplot.plot(min_gamma + numpy.argmax(ref_sa_plus_petfold_mccs), max(ref_sa_plus_petfold_mccs), label = "PETfold", marker = "^", linestyle = "-.", markerfacecolor = white, markeredgecolor = color_palette[2])
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("Matthews correlation coefficient")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_mccs_on_css_estimation_ref_sa.eps", bbox_inches = "tight")
  pyplot.clf()
  # Figure for LocARNA.
  line_1, = pyplot.plot(mafft_xinsi_plus_consalifold_ppvs, mafft_xinsi_plus_consalifold_senss, label = "MAFFT X-INS-i + ConsAlifold", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(ref_sa_plus_consalifold_ppvs, ref_sa_plus_consalifold_senss, label = "Reference + ConsAlifold", marker = "s", linestyle = "-")
  line_3, = pyplot.plot(centroidhomfold_ppvs, centroidhomfold_senss, label = "CentroidHomfold", marker = "^", linestyle = "--")
  line_4, = pyplot.plot(locarna_ppv, locarna_sens, label = "LocARNA", marker = "v", linestyle = "-")
  line_5, = pyplot.plot(raf_ppv, raf_sens, label = "RAF", marker = "d", linestyle = "-")
  line_6, = pyplot.plot(turbofold_ppvs, turbofold_senss, label = "TurboFold", marker = "p", linestyle = "-.")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4, line_5, line_6], loc = "lower left")
  pyplot.xlabel("Precision")
  pyplot.ylabel("Recall")
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/pr_curves_on_ss_estimation_other_types.eps", bbox_inches = "tight")
  pyplot.clf()
  pyplot.figure()
  line_1, = pyplot.plot(mafft_xinsi_plus_consalifold_fprs, mafft_xinsi_plus_consalifold_senss, label = "MAFFT X-INS-i + ConsAlifold", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(ref_sa_plus_consalifold_fprs, ref_sa_plus_consalifold_senss, label = "Reference + ConsAlifold", marker = "s", linestyle = "-")
  line_3, = pyplot.plot(centroidhomfold_fprs, centroidhomfold_senss, label = "CentroidHomfold", marker = "^", linestyle = "--")
  line_4, = pyplot.plot(locarna_fpr, locarna_sens, label = "LocARNA", marker = "v", linestyle = "-")
  line_5, = pyplot.plot(raf_fpr, raf_sens, label = "RAF", marker = "d", linestyle = "-")
  line_6, = pyplot.plot(turbofold_fprs, turbofold_senss, label = "TurboFold", marker = "p", linestyle = "-.")
  pyplot.xlabel("Fall-out")
  pyplot.ylabel("Recall")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/roc_curves_on_ss_estimation_other_types.eps", bbox_inches = "tight")
  pyplot.clf()
  gammas = [i for i in range(min_gamma, max_gamma + 1)]
  pyplot.figure()
  line_1, = pyplot.plot(gammas, mafft_xinsi_plus_consalifold_f1_scores, label = "MAFFT X-INS-i + ConsAlifold", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, ref_sa_plus_consalifold_f1_scores, label = "Reference + ConsAlifold", marker = "s", linestyle = "-")
  line_3, = pyplot.plot(gammas, centroidhomfold_f1_scores, label = "CentroidHomfold", marker = "^", linestyle = "--")
  line_4, = pyplot.plot(-3., locarna_f1_score, label = "LocARNA", marker = "v", linestyle = "-", zorder = 9)
  line_5, = pyplot.plot(-3., raf_f1_score, label = "RAF", marker = "d", linestyle = "-", zorder = 10)
  line_6, = pyplot.plot(gammas, turbofold_f1_scores, label = "TurboFold", marker = "p", linestyle = "-.")
  line_7, = pyplot.plot(min_gamma + numpy.argmax(mafft_xinsi_plus_consalifold_f1_scores), max(mafft_xinsi_plus_consalifold_f1_scores), label = "MAFFT X-INS-i + ConsAlifold", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[0])
  line_8, = pyplot.plot(min_gamma + numpy.argmax(ref_sa_plus_consalifold_f1_scores), max(ref_sa_plus_consalifold_f1_scores), label = "Reference + ConsAlifold", marker = "s", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[1])
  line_9, = pyplot.plot(min_gamma + numpy.argmax(centroidhomfold_f1_scores), max(centroidhomfold_f1_scores), label = "CentroidHomfold", marker = "^", linestyle = "--", markerfacecolor = white, markeredgecolor = color_palette[2])
  line_10, = pyplot.plot(min_gamma + numpy.argmax(turbofold_f1_scores), max(turbofold_f1_scores), label = "TurboFold", marker = "p", linestyle = "-.", markerfacecolor = white, markeredgecolor = color_palette[5])
  pyplot.legend(handles = [line_7, line_8, line_9, line_10], loc = "lower right")
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("F1 score")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_f1_scores_on_ss_estimation_other_types.eps", bbox_inches = "tight")
  pyplot.clf()
  pyplot.figure()
  line_1, = pyplot.plot(gammas, mafft_xinsi_plus_consalifold_mccs, label = "MAFFT X-INS-i + ConsAlifold", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(gammas, ref_sa_plus_consalifold_mccs, label = "Reference + ConsAlifold", marker = "s", linestyle = "-")
  line_3, = pyplot.plot(gammas, centroidhomfold_mccs, label = "CentroidHomfold", marker = "^", linestyle = "--")
  line_4, = pyplot.plot(-3., locarna_mcc, label = "LocARNA", marker = "v", linestyle = "-", zorder = 9)
  line_5, = pyplot.plot(-3., raf_mcc, label = "RAF", marker = "d", linestyle = "-", zorder = 10)
  line_6, = pyplot.plot(gammas, turbofold_mccs, label = "TurboFold", marker = "p", linestyle = "-.")
  line_7, = pyplot.plot(min_gamma + numpy.argmax(mafft_xinsi_plus_consalifold_mccs), max(mafft_xinsi_plus_consalifold_mccs), label = "MAFFT X-INS-i + ConsAlifold", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[0])
  line_8, = pyplot.plot(min_gamma + numpy.argmax(ref_sa_plus_consalifold_mccs), max(ref_sa_plus_consalifold_mccs), label = "Reference + ConsAlifold", marker = "s", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[1])
  line_9, = pyplot.plot(min_gamma + numpy.argmax(centroidhomfold_mccs), max(centroidhomfold_mccs), label = "CentroidHomfold", marker = "^", linestyle = "--", markerfacecolor = white, markeredgecolor = color_palette[2])
  line_10, = pyplot.plot(min_gamma + numpy.argmax(turbofold_mccs), max(turbofold_mccs), label = "TurboFold", marker = "p", linestyle = "-.", markerfacecolor = white, markeredgecolor = color_palette[5])
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("Matthews correlation coefficient")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_mccs_on_ss_estimation_other_types.eps", bbox_inches = "tight")
  pyplot.clf()

def get_metrics(bin_counts):
  (tp, tn, fp, fn) = bin_counts
  ppv = get_ppv(tp, fp)
  sens = get_sens(tp, fn)
  fpr = get_fpr(tn, fp)
  f1_score = get_f1_score(ppv, sens)
  mcc = get_mcc(tp, tn, fp, fn)
  return ppv, sens, fpr, f1_score, mcc

def get_bin_counts(params):
  rna_seq_lens, estimated_css, ref_css = params
  num_of_rnas = len(rna_seq_lens)
  tp = fp = tn = fn = 0
  for m in range(0, num_of_rnas):
    sub_estimated_css = estimated_css[m]
    sub_ref_css = ref_css[m]
    rna_seq_len_1 = rna_seq_lens[m]
    for i in range(0, rna_seq_len_1):
      for j in range(i + 1, rna_seq_len_1):
        estimated_bin = (i, j) in sub_estimated_css
        ref_bin = (i, j) in sub_ref_css
        if estimated_bin == ref_bin:
          if estimated_bin == True:
            tp += 1
          else:
            tn += 1
        else:
          if estimated_bin == True:
            fp += 1
          else:
            fn += 1
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

def get_sss(ss_file_path):
  sss = []
  ss_strings = []
  ss_strings = [rec.seq for rec in SeqIO.parse(ss_file_path, "fasta")]
  sss = []
  for (i, ss_string) in enumerate(ss_strings):
    sss.append({})
    for (left, right) in bracket_pairs:
      stack = []
      for (j, char) in enumerate(ss_string):
        if char == left:
          stack.append(j)
        elif char == right:
          pos = stack.pop()
          sss[i][(pos, j)] = True
  return sss

if __name__ == "__main__":
  main()
