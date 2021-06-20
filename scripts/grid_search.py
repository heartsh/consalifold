#! /usr/bin/env python

import utils
from Bio import SeqIO
import numpy
import seaborn
from matplotlib import pyplot
import os
import multiprocessing
import time
import datetime
import shutil
from Bio import AlignIO
import math
from math import sqrt

taus = numpy.arange(0, 1.1, 0.1).tolist()
seaborn.set()
pyplot.rcParams['legend.handlelength'] = 0
pyplot.rcParams['legend.fontsize'] = "large"
color_palette = seaborn.color_palette()
min_gamma = -4
max_gamma = 10
white = "#F2F2F2"

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  num_of_threads = multiprocessing.cpu_count()
  consalifold_params = []
  rna_seq_dir_path = asset_dir_path + "/compiled_rna_fams_valid"
  ref_sa_dir_path = asset_dir_path + "/ref_sas_valid"
  consalifold_dir_path = asset_dir_path + "/grid_search"
  if not os.path.isdir(consalifold_dir_path):
    os.mkdir(consalifold_dir_path)
  sub_thread_num = 4
  for tau in taus:
    tau_str = "%.1f" % tau
    consalifold_sub_dir_path = os.path.join(consalifold_dir_path, tau_str)
    if not os.path.isdir(consalifold_sub_dir_path):
      os.mkdir(consalifold_sub_dir_path)
    for rna_seq_file in os.listdir(rna_seq_dir_path):
      if not rna_seq_file.endswith(".fa"):
        continue
      rna_seq_file_path = os.path.join(rna_seq_dir_path, rna_seq_file)
      (rna_family_name, extension) = os.path.splitext(rna_seq_file)
      ref_sa_file_path = os.path.join(ref_sa_dir_path, rna_family_name + ".aln")
      consalifold_output_dir_path = os.path.join(consalifold_sub_dir_path, rna_family_name)
      if not os.path.isdir(consalifold_output_dir_path):
        os.mkdir(consalifold_output_dir_path)
      consalifold_command = "consalifold -t " + str(sub_thread_num) + " -i " + ref_sa_file_path + " -o " + consalifold_output_dir_path + " --mix_weight " + str(tau)
      consalifold_params.insert(0, consalifold_command)
  # ConsAliFold's execution.
  pool = multiprocessing.Pool(int(num_of_threads / sub_thread_num))
  pool.map(utils.run_command, consalifold_params)
  max_consalifold_mccs = [0] * len(taus)
  gammas = [2. ** i for i in range(min_gamma, max_gamma + 1)]
  for i in range(len(taus)):
    tau = taus[i]
    tau_str = "%.1f" % tau
    consalifold_sub_dir_path = os.path.join(consalifold_dir_path, tau_str)
    for gamma in gammas:
      consalifold_count_params = []
      gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
      for rna_fam_file in os.listdir(rna_seq_dir_path):
        if not rna_fam_file.endswith(".fa"):
          continue
        rna_seq_file_path = os.path.join(rna_seq_dir_path, rna_fam_file)
        rna_seq_lens = [len(rna_seq.seq) for rna_seq in SeqIO.parse(rna_seq_file_path, "fasta")]
        (rna_fam_name, extension) = os.path.splitext(rna_fam_file)
        ref_css_file_path = os.path.join(ref_sa_dir_path, rna_fam_name + ".sth")
        ref_css = utils.get_css(ref_css_file_path)
        consalifold_estimated_css_dir_path = os.path.join(consalifold_sub_dir_path, rna_fam_name)
        consalifold_estimated_css_file_path = os.path.join(consalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
        estimated_css = utils.get_css(consalifold_estimated_css_file_path)
        consalifold_count_params.insert(0, (rna_seq_lens, estimated_css, ref_css))
      results = pool.map(get_bin_counts, consalifold_count_params)
      _, _, _, _, mcc = get_metrics(final_sum(results))
      max_consalifold_mcc = max_consalifold_mccs[i]
      if mcc > max_consalifold_mcc:
        max_consalifold_mccs[i] = mcc
  best_tau = taus[numpy.argmax(max_consalifold_mccs)]
  print("The best mixing coefficient based on maximum MCCs = %.1f." % best_tau)
  line_1, = pyplot.plot(taus, max_consalifold_mccs, marker = "o", linestyle = "-")
  line_2, = pyplot.plot(best_tau, max(max_consalifold_mccs), label = "Best probability mixing coefficient", marker = "o", linestyle = "-", markerfacecolor = white, markeredgecolor = color_palette[0])
  pyplot.legend(handles = [line_2], loc = "lower left")
  pyplot.xlabel("Posterior probability mixing coefficient")
  pyplot.ylabel("Maximum Matthews correlation coefficient")
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/mixing_coefficient_grid_search.eps", bbox_inches = "tight")
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

if __name__ == "__main__":
  main()
