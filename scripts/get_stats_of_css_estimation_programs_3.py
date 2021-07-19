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
  locarna_plus_consalifold_ppvs = []
  locarna_plus_consalifold_senss = []
  locarna_plus_consalifold_fprs = []
  locarna_plus_consalifold_f1_scores = []
  locarna_plus_consalifold_mccs = []
  raf_plus_consalifold_ppvs = []
  raf_plus_consalifold_senss = []
  raf_plus_consalifold_fprs = []
  raf_plus_consalifold_f1_scores = []
  raf_plus_consalifold_mccs = []
  locarna_ppv = locarna_sens = locarna_fpr = locarna_f1_score = locarna_mcc = 0.
  raf_ppv = raf_sens = raf_fpr = raf_f1_score = raf_mcc = 0.
  gammas = [2. ** i for i in range(min_gamma, max_gamma + 1)]
  rna_fam_dir_path = asset_dir_path + "/compiled_rna_fams_test"
  ref_sa_dir_path = asset_dir_path + "/ref_sas_test"
  locarna_css_dir_path = asset_dir_path + "/locarna"
  raf_css_dir_path = asset_dir_path + "/raf"
  locarna_plus_consalifold_css_dir_path = asset_dir_path + "/locarna_plus_consalifold"
  raf_plus_consalifold_css_dir_path = asset_dir_path + "/raf_plus_consalifold"
  pool = multiprocessing.Pool(num_of_threads)
  for gamma in gammas:
    locarna_count_params = []
    raf_count_params = []
    locarna_plus_consalifold_count_params = []
    raf_plus_consalifold_count_params = []
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
      locarna_plus_consalifold_estimated_css_dir_path = os.path.join(locarna_plus_consalifold_css_dir_path, rna_fam_name)
      raf_plus_consalifold_estimated_css_dir_path = os.path.join(raf_plus_consalifold_css_dir_path, rna_fam_name)
      locarna_plus_consalifold_estimated_css_file_path = os.path.join(locarna_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(locarna_plus_consalifold_estimated_css_file_path)
      locarna_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      raf_plus_consalifold_estimated_css_file_path = os.path.join(raf_plus_consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css = utils.get_css(raf_plus_consalifold_estimated_css_file_path)
      raf_plus_consalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      if gamma == 1.:
        locarna_estimated_css_file_path = os.path.join(locarna_css_dir_path, rna_fam_name + "/results/result.stk")
        estimated_css = utils.get_css(locarna_estimated_css_file_path)
        locarna_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
        raf_estimated_css_file_path = os.path.join(raf_css_dir_path, rna_fam_name + ".sth")
        estimated_css = utils.get_css(raf_estimated_css_file_path)
        raf_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
    results = pool.map(get_bin_counts, locarna_plus_consalifold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    locarna_plus_consalifold_ppvs.insert(0, ppv)
    locarna_plus_consalifold_senss.insert(0, sens)
    locarna_plus_consalifold_fprs.insert(0, fpr)
    locarna_plus_consalifold_f1_scores.append(f1_score)
    locarna_plus_consalifold_mccs.append(mcc)
    results = pool.map(get_bin_counts, raf_plus_consalifold_count_params)
    ppv, sens, fpr, f1_score, mcc = get_metrics(final_sum(results))
    raf_plus_consalifold_ppvs.insert(0, ppv)
    raf_plus_consalifold_senss.insert(0, sens)
    raf_plus_consalifold_fprs.insert(0, fpr)
    raf_plus_consalifold_f1_scores.append(f1_score)
    raf_plus_consalifold_mccs.append(mcc)
    if gamma == 1.:
      results = pool.map(get_bin_counts, locarna_count_params)
      locarna_ppv, locarna_sens, locarna_fpr, locarna_f1_score, locarna_mcc = get_metrics(final_sum(results))
      results = pool.map(get_bin_counts, raf_count_params)
      raf_ppv, raf_sens, raf_fpr, raf_f1_score, raf_mcc = get_metrics(final_sum(results))
  line_1, = pyplot.plot(locarna_ppv, locarna_sens, label = "LocARNA", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(raf_ppv, raf_sens, label = "RAF", marker = "v", linestyle = "--")
  line_3, = pyplot.plot(locarna_plus_consalifold_ppvs, locarna_plus_consalifold_senss, label = "LocARNA + ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_4, = pyplot.plot(raf_plus_consalifold_ppvs, raf_plus_consalifold_senss, label = "RAF + ConsAlifold (ConsProb)", marker = "v", linestyle = "--")
  pyplot.xlabel("Precision")
  pyplot.ylabel("Recall")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4], loc = "lower left")
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/pr_curves_on_css_estimation_struct_aligners.eps", bbox_inches = "tight")
  pyplot.clf()
  pyplot.figure()
  line_1, = pyplot.plot(locarna_fpr, locarna_sens, label = "LocARNA", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(raf_fpr, raf_sens, label = "RAF", marker = "v", linestyle = "--")
  line_3, = pyplot.plot(locarna_plus_consalifold_fprs, locarna_plus_consalifold_senss, label = "LocARNA + ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_4, = pyplot.plot(raf_plus_consalifold_fprs, raf_plus_consalifold_senss, label = "RAF + ConsAlifold (ConsProb)", marker = "v", linestyle = "--")
  pyplot.xlabel("Fall-out")
  pyplot.ylabel("Recall")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/roc_curves_on_css_estimation_struct_aligners.eps", bbox_inches = "tight")
  pyplot.clf()
  gammas = [i for i in range(min_gamma, max_gamma + 1)]
  pyplot.figure()
  line_1, = pyplot.plot(-3., locarna_f1_score, label = "LocARNA", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(-3., raf_f1_score, label = "RAF", marker = "v", linestyle = "--")
  line_3, = pyplot.plot(gammas, locarna_plus_consalifold_f1_scores, label = "LocARNA + ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_4, = pyplot.plot(gammas, raf_plus_consalifold_f1_scores, label = "RAF + ConsAlifold (ConsProb)", marker = "v", linestyle = "--")
  line_5, = pyplot.plot(min_gamma + numpy.argmax(locarna_plus_consalifold_f1_scores), max(locarna_plus_consalifold_f1_scores), label = "LocARNA + ConsAlifold (ConsProb)", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[2])
  line_6, = pyplot.plot(min_gamma + numpy.argmax(raf_plus_consalifold_f1_scores), max(raf_plus_consalifold_f1_scores), label = "RAF + ConsAlifold (ConsProb)", marker = "v", linestyle = "--", markerfacecolor = white, markeredgecolor = color_palette[3])
  pyplot.legend(handles = [line_5, line_6], loc = "lower right")
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("F1 score")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_f1_scores_on_css_estimation_struct_aligners.eps", bbox_inches = "tight")
  pyplot.clf()
  pyplot.figure()
  line_1, = pyplot.plot(-3., locarna_mcc, label = "LocARNA", marker = "o", linestyle = "-")
  line_2, = pyplot.plot(-3., raf_mcc, label = "RAF", marker = "v", linestyle = "--")
  line_3, = pyplot.plot(gammas, locarna_plus_consalifold_mccs, label = "LocARNA + ConsAlifold (ConsProb)", marker = "o", linestyle = "-")
  line_4, = pyplot.plot(gammas, raf_plus_consalifold_mccs, label = "RAF + ConsAlifold (ConsProb)", marker = "v", linestyle = "--")
  line_5, = pyplot.plot(min_gamma + numpy.argmax(locarna_plus_consalifold_mccs), max(locarna_plus_consalifold_mccs), label = "LocARNA + ConsAlifold (ConsProb)", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[2])
  line_6, = pyplot.plot(min_gamma + numpy.argmax(raf_plus_consalifold_mccs), max(raf_plus_consalifold_mccs), label = "RAF + ConsAlifold (ConsProb)", marker = "v", linestyle = "--", markerfacecolor = white, markeredgecolor = color_palette[3])
  pyplot.xlabel("$\log_2 \gamma$")
  pyplot.ylabel("Matthews correlation coefficient")
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/gammas_vs_mccs_on_css_estimation_struct_aligners.eps", bbox_inches = "tight")
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
