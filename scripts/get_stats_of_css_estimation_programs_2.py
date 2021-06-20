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
import pandas
from scipy import stats

seaborn.set()
# pyplot.rcParams['legend.handlelength'] = 0
# pyplot.rcParams['legend.fontsize'] = "x-large"
# pyplot.rcParams['legend.fontsize'] = "large"
# color_palette = seaborn.color_palette()
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
  posterior_probcons_plus_consalifold_ppvs = []
  posterior_probcons_plus_consalifold_senss = []
  posterior_probcons_plus_consalifold_fprs = []
  posterior_probcons_plus_consalifold_f1_scores = []
  posterior_probcons_plus_consalifold_mccs = []
  posterior_clustalw_plus_consalifold_ppvs = []
  posterior_clustalw_plus_consalifold_senss = []
  posterior_clustalw_plus_consalifold_fprs = []
  posterior_clustalw_plus_consalifold_f1_scores = []
  posterior_clustalw_plus_consalifold_mccs = []
  posterior_mafft_plus_consalifold_ppvs = []
  posterior_mafft_plus_consalifold_senss = []
  posterior_mafft_plus_consalifold_fprs = []
  posterior_mafft_plus_consalifold_f1_scores = []
  posterior_mafft_plus_consalifold_mccs = []
  posterior_mafft_xinsi_plus_consalifold_ppvs = []
  posterior_mafft_xinsi_plus_consalifold_senss = []
  posterior_mafft_xinsi_plus_consalifold_fprs = []
  posterior_mafft_xinsi_plus_consalifold_f1_scores = []
  posterior_mafft_xinsi_plus_consalifold_mccs = []
  posterior_ref_sa_plus_consalifold_ppvs = []
  posterior_ref_sa_plus_consalifold_senss = []
  posterior_ref_sa_plus_consalifold_fprs = []
  posterior_ref_sa_plus_consalifold_f1_scores = []
  posterior_ref_sa_plus_consalifold_mccs = []
  gammas = [2. ** i for i in range(min_gamma, max_gamma + 1)]
  rna_fam_dir_path = asset_dir_path + "/compiled_rna_fams_test"
  ref_sa_dir_path = asset_dir_path + "/ref_sas_test"
  mafft_plus_consalifold_css_dir_path = asset_dir_path + "/mafft_plus_consalifold"
  probcons_plus_consalifold_css_dir_path = asset_dir_path + "/probcons_plus_consalifold"
  clustalw_plus_consalifold_css_dir_path = asset_dir_path + "/clustalw_plus_consalifold"
  mafft_xinsi_plus_consalifold_css_dir_path = asset_dir_path + "/mafft_xinsi_plus_consalifold"
  ref_sa_plus_consalifold_css_dir_path = asset_dir_path + "/ref_sa_plus_consalifold"
  posterior_probcons_plus_consalifold_css_dir_path = asset_dir_path + "/posterior_probcons_plus_consalifold"
  posterior_clustalw_plus_consalifold_css_dir_path = asset_dir_path + "/posterior_clustalw_plus_consalifold"
  posterior_mafft_plus_consalifold_css_dir_path = asset_dir_path + "/posterior_mafft_plus_consalifold"
  posterior_mafft_xinsi_plus_consalifold_css_dir_path = asset_dir_path + "/posterior_mafft_xinsi_plus_consalifold"
  posterior_ref_sa_plus_consalifold_css_dir_path = asset_dir_path + "/posterior_ref_sa_plus_consalifold"
  pool = multiprocessing.Pool(num_of_threads)
  for gamma in gammas:
    mafft_plus_consalifold_count_params = []
    clustalw_plus_consalifold_count_params = []
    mafft_xinsi_plus_consalifold_count_params = []
    ref_sa_plus_consalifold_count_params = []
    probcons_plus_consalifold_count_params = []
    posterior_probcons_plus_consalifold_count_params = []
    posterior_clustalw_plus_consalifold_count_params = []
    posterior_mafft_plus_consalifold_count_params = []
    posterior_mafft_xinsi_plus_consalifold_count_params = []
    posterior_ref_sa_plus_consalifold_count_params = []
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
      posterior_probcons_plus_consalifold_estimated_css_dir_path = os.path.join(posterior_probcons_plus_consalifold_css_dir_path, rna_fam_name)
      posterior_clustalw_plus_consalifold_estimated_css_dir_path = os.path.join(posterior_clustalw_plus_consalifold_css_dir_path, rna_fam_name)
      posterior_mafft_plus_consalifold_estimated_css_dir_path = os.path.join(posterior_mafft_plus_consalifold_css_dir_path, rna_fam_name)
      posterior_mafft_xinsi_plus_consalifold_estimated_css_dir_path = os.path.join(posterior_mafft_xinsi_plus_consalifold_css_dir_path, rna_fam_name)
      posterior_ref_sa_plus_consalifold_estimated_css_dir_path = os.path.join(posterior_ref_sa_plus_consalifold_css_dir_path, rna_fam_name)
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
      posterior_probcons_plus_consalifold_estimated_css_file_path = os.path.join(posterior_probcons_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(posterior_probcons_plus_consalifold_estimated_css_file_path)
      posterior_probcons_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      posterior_clustalw_plus_consalifold_estimated_css_file_path = os.path.join(posterior_clustalw_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(posterior_clustalw_plus_consalifold_estimated_css_file_path)
      posterior_clustalw_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      posterior_mafft_plus_consalifold_estimated_css_file_path = os.path.join(posterior_mafft_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(posterior_mafft_plus_consalifold_estimated_css_file_path)
      posterior_mafft_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      posterior_mafft_xinsi_plus_consalifold_estimated_css_file_path = os.path.join(posterior_mafft_xinsi_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(posterior_mafft_xinsi_plus_consalifold_estimated_css_file_path)
      posterior_mafft_xinsi_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      posterior_ref_sa_plus_consalifold_estimated_css_file_path = os.path.join(posterior_ref_sa_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(posterior_ref_sa_plus_consalifold_estimated_css_file_path)
      posterior_ref_sa_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
    results = pool.map(get_bin_counts, probcons_plus_consalifold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    probcons_plus_consalifold_ppvs.insert(0, ppv)
    probcons_plus_consalifold_senss.insert(0, sens)
    probcons_plus_consalifold_fprs.insert(0, fpr)
    probcons_plus_consalifold_f1_scores.append(f1_score)
    probcons_plus_consalifold_mccs.append(mcc)
    if True:
      results = pool.map(get_bin_counts, clustalw_plus_consalifold_count_params)
      ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
      clustalw_plus_consalifold_ppvs.insert(0, ppv)
      clustalw_plus_consalifold_senss.insert(0, sens)
      clustalw_plus_consalifold_fprs.insert(0, fpr)
      clustalw_plus_consalifold_f1_scores.append(f1_score)
      clustalw_plus_consalifold_mccs.append(mcc)
      results = pool.map(get_bin_counts, mafft_plus_consalifold_count_params)
      ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
      mafft_plus_consalifold_ppvs.insert(0, ppv)
      mafft_plus_consalifold_senss.insert(0, sens)
      mafft_plus_consalifold_fprs.insert(0, fpr)
      mafft_plus_consalifold_f1_scores.append(f1_score)
      mafft_plus_consalifold_mccs.append(mcc)
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
    results = pool.map(get_bin_counts, posterior_probcons_plus_consalifold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    posterior_probcons_plus_consalifold_ppvs.insert(0, ppv)
    posterior_probcons_plus_consalifold_senss.insert(0, sens)
    posterior_probcons_plus_consalifold_fprs.insert(0, fpr)
    posterior_probcons_plus_consalifold_f1_scores.append(f1_score)
    posterior_probcons_plus_consalifold_mccs.append(mcc)
    if True:
      results = pool.map(get_bin_counts, posterior_clustalw_plus_consalifold_count_params)
      ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
      posterior_clustalw_plus_consalifold_ppvs.insert(0, ppv)
      posterior_clustalw_plus_consalifold_senss.insert(0, sens)
      posterior_clustalw_plus_consalifold_fprs.insert(0, fpr)
      posterior_clustalw_plus_consalifold_f1_scores.append(f1_score)
      posterior_clustalw_plus_consalifold_mccs.append(mcc)
      results = pool.map(get_bin_counts, posterior_mafft_plus_consalifold_count_params)
      ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
      posterior_mafft_plus_consalifold_ppvs.insert(0, ppv)
      posterior_mafft_plus_consalifold_senss.insert(0, sens)
      posterior_mafft_plus_consalifold_fprs.insert(0, fpr)
      posterior_mafft_plus_consalifold_f1_scores.append(f1_score)
      posterior_mafft_plus_consalifold_mccs.append(mcc)
      results = pool.map(get_bin_counts, posterior_mafft_xinsi_plus_consalifold_count_params)
      ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
      posterior_mafft_xinsi_plus_consalifold_ppvs.insert(0, ppv)
      posterior_mafft_xinsi_plus_consalifold_senss.insert(0, sens)
      posterior_mafft_xinsi_plus_consalifold_fprs.insert(0, fpr)
      posterior_mafft_xinsi_plus_consalifold_f1_scores.append(f1_score)
      posterior_mafft_xinsi_plus_consalifold_mccs.append(mcc)
      results = pool.map(get_bin_counts, posterior_ref_sa_plus_consalifold_count_params)
      ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
      posterior_ref_sa_plus_consalifold_ppvs.insert(0, ppv)
      posterior_ref_sa_plus_consalifold_senss.insert(0, sens)
      posterior_ref_sa_plus_consalifold_fprs.insert(0, fpr)
      posterior_ref_sa_plus_consalifold_f1_scores.append(f1_score)
      posterior_ref_sa_plus_consalifold_mccs.append(mcc)
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  consalifold_avg_mccs = [numpy.mean(clustalw_plus_consalifold_mccs), numpy.mean(mafft_plus_consalifold_mccs), numpy.mean(probcons_plus_consalifold_mccs), numpy.mean(mafft_xinsi_plus_consalifold_mccs), numpy.mean(ref_sa_plus_consalifold_mccs)]
  posterior_consalifold_avg_mccs = [numpy.mean(posterior_clustalw_plus_consalifold_mccs), numpy.mean(posterior_mafft_plus_consalifold_mccs), numpy.mean(posterior_probcons_plus_consalifold_mccs), numpy.mean(posterior_mafft_xinsi_plus_consalifold_mccs), numpy.mean(posterior_ref_sa_plus_consalifold_mccs)]
  avg_mccs = consalifold_avg_mccs + posterior_consalifold_avg_mccs
  data = {"Average Matthews correlation coefficient": avg_mccs, "Alignment probability inference method": ["ConsProb"] * 5 + ["LocARNA-P + our PCT"] * 5, "Sequence alignment source": ["ClustalW", "MAFFT", "ProbCons-RNA ", "MAFFT X-INS-i", "Reference"] * 2}
  data_frame = pandas.DataFrame(data = data)
  ax = seaborn.barplot(x = "Sequence alignment source", y = "Average Matthews correlation coefficient", hue = "Alignment probability inference method", data = data_frame)
  # pyplot.ylim(0, 0.8)
  ax.legend_.remove()
  fig = ax.get_figure()
  fig.tight_layout()
  fig.savefig(image_dir_path + "/consalifold_model_comparison_mcc.eps", bbox_inches = "tight")
  fig.clf()
  consalifold_avg_f1_scores = [numpy.mean(clustalw_plus_consalifold_f1_scores), numpy.mean(mafft_plus_consalifold_f1_scores), numpy.mean(probcons_plus_consalifold_f1_scores), numpy.mean(mafft_xinsi_plus_consalifold_f1_scores), numpy.mean(ref_sa_plus_consalifold_f1_scores)]
  posterior_consalifold_avg_f1_scores = [numpy.mean(posterior_clustalw_plus_consalifold_f1_scores), numpy.mean(posterior_mafft_plus_consalifold_f1_scores), numpy.mean(posterior_probcons_plus_consalifold_f1_scores), numpy.mean(posterior_mafft_xinsi_plus_consalifold_f1_scores), numpy.mean(posterior_ref_sa_plus_consalifold_f1_scores)]
  avg_f1_scores = consalifold_avg_f1_scores + posterior_consalifold_avg_f1_scores
  data = {"Average F1 score": avg_f1_scores, "Alignment probability inference method": ["ConsProb"] * 5 + ["LocARNA-P + our PCT"] * 5, "Sequence alignment source": ["ClustalW", "MAFFT", "ProbCons-RNA ", "MAFFT X-INS-i", "Reference"] * 2}
  data_frame = pandas.DataFrame(data = data)
  ax = seaborn.barplot(x = "Sequence alignment source", y = "Average F1 score", hue = "Alignment probability inference method", data = data_frame)
  pyplot.ylim(0, 0.75)
  ax.legend(loc = "upper left")
  fig = ax.get_figure()
  fig.tight_layout()
  fig.savefig(image_dir_path + "/consalifold_model_comparison_f1_score.eps", bbox_inches = "tight")
  fig.clf()
  consalifold_mccs = clustalw_plus_consalifold_mccs + mafft_plus_consalifold_mccs + probcons_plus_consalifold_mccs + mafft_xinsi_plus_consalifold_mccs + ref_sa_plus_consalifold_mccs
  posterior_consalifold_mccs = posterior_clustalw_plus_consalifold_mccs + posterior_mafft_plus_consalifold_mccs + posterior_probcons_plus_consalifold_mccs + posterior_mafft_xinsi_plus_consalifold_mccs + posterior_ref_sa_plus_consalifold_mccs
  consalifold_f1_scores = clustalw_plus_consalifold_f1_scores + mafft_plus_consalifold_f1_scores + probcons_plus_consalifold_f1_scores + mafft_xinsi_plus_consalifold_f1_scores + ref_sa_plus_consalifold_f1_scores
  posterior_consalifold_f1_scores = posterior_clustalw_plus_consalifold_f1_scores + posterior_mafft_plus_consalifold_f1_scores + posterior_probcons_plus_consalifold_f1_scores + posterior_mafft_xinsi_plus_consalifold_f1_scores + posterior_ref_sa_plus_consalifold_f1_scores
  print("ClustalW (MCCs used):", stats.ttest_rel(clustalw_plus_consalifold_mccs, posterior_clustalw_plus_consalifold_mccs))
  print("MAFFT (MCCs used):", stats.ttest_rel(mafft_plus_consalifold_mccs, posterior_mafft_plus_consalifold_mccs))
  print("ProbCons (MCCs used):", stats.ttest_rel(probcons_plus_consalifold_mccs, posterior_probcons_plus_consalifold_mccs))
  print("MAFFT X-INS-i (MCCs used):", stats.ttest_rel(mafft_xinsi_plus_consalifold_mccs, posterior_mafft_xinsi_plus_consalifold_mccs))
  print("Reference (MCCs used):", stats.ttest_rel(ref_sa_plus_consalifold_mccs, posterior_ref_sa_plus_consalifold_mccs))
  print("Overall (MCCs used):", stats.ttest_rel(consalifold_mccs, posterior_consalifold_mccs))
  print("ClustalW (F1 scores used):", stats.ttest_rel(clustalw_plus_consalifold_f1_scores, posterior_clustalw_plus_consalifold_f1_scores))
  print("MAFFT (F1 scores used):", stats.ttest_rel(mafft_plus_consalifold_f1_scores, posterior_mafft_plus_consalifold_f1_scores))
  print("ProbCons (F1 scores used):", stats.ttest_rel(probcons_plus_consalifold_f1_scores, posterior_probcons_plus_consalifold_f1_scores))
  print("MAFFT X-INS-i (F1 scores used):", stats.ttest_rel(mafft_xinsi_plus_consalifold_f1_scores, posterior_mafft_xinsi_plus_consalifold_f1_scores))
  print("Reference (F1 scores used):", stats.ttest_rel(ref_sa_plus_consalifold_f1_scores, posterior_ref_sa_plus_consalifold_f1_scores))
  print("Overall (F1 scores used):", stats.ttest_rel(consalifold_f1_scores, posterior_consalifold_f1_scores))

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
