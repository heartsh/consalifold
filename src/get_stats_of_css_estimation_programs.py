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
  mafft_ginsi_plus_neoalifold_css_dir_path = asset_dir_path + "/mafft_ginsi_plus_neoalifold"
  mafft_xinsi_plus_neoalifold_css_dir_path = asset_dir_path + "/mafft_xinsi_plus_neoalifold"
  ref_sa_plus_neoalifold_css_dir_path = asset_dir_path + "/ref_sa_plus_neoalifold"
  mafft_ginsi_plus_centroidalifold_css_dir_path = asset_dir_path + "/mafft_ginsi_plus_centroidalifold"
  mafft_xinsi_plus_centroidalifold_css_dir_path = asset_dir_path + "/mafft_xinsi_plus_centroidalifold"
  ref_sa_plus_centroidalifold_css_dir_path = asset_dir_path + "/ref_sa_plus_centroidalifold"
  rna_fam_dir_path = asset_dir_path + "/sampled_rna_fams"
  mafft_ginsi_plus_neoalifold_ppvs = []
  mafft_ginsi_plus_neoalifold_senss = []
  mafft_ginsi_plus_neoalifold_fprs = []
  mafft_xinsi_plus_neoalifold_ppvs = []
  mafft_xinsi_plus_neoalifold_senss = []
  mafft_xinsi_plus_neoalifold_fprs = []
  ref_sa_plus_neoalifold_ppvs = []
  ref_sa_plus_neoalifold_senss = []
  ref_sa_plus_neoalifold_fprs = []
  # mafft_ginsi_plus_centroidalifold_ppvs = []
  # mafft_ginsi_plus_centroidalifold_senss = []
  # mafft_ginsi_plus_centroidalifold_fprs = []
  # mafft_xinsi_plus_centroidalifold_ppvs = []
  # mafft_xinsi_plus_centroidalifold_senss = []
  # mafft_xinsi_plus_centroidalifold_fprs = []
  # ref_sa_plus_centroidalifold_ppvs = []
  # ref_sa_plus_centroidalifold_senss = []
  # ref_sa_plus_centroidalifold_fprs = []
  gammas = [2. ** i for i in range(-7, 11)]
  for gamma in gammas:
    gamma_str = str(gamma)
    mafft_ginsi_plus_neoalifold_tp = mafft_ginsi_plus_neoalifold_tn = mafft_ginsi_plus_neoalifold_fp = mafft_ginsi_plus_neoalifold_fn = 0.
    mafft_xinsi_plus_neoalifold_tp = mafft_xinsi_plus_neoalifold_tn = mafft_xinsi_plus_neoalifold_fp = mafft_xinsi_plus_neoalifold_fn = 0.
    ref_sa_plus_neoalifold_tp = ref_sa_plus_neoalifold_tn = ref_sa_plus_neoalifold_fp = ref_sa_plus_neoalifold_fn = 0.
    # mafft_ginsi_plus_centroidalifold_tp = mafft_ginsi_plus_centroidalifold_tn = mafft_ginsi_plus_centroidalifold_fp = mafft_ginsi_plus_centroidalifold_fn = 0.
    # mafft_xinsi_plus_centroidalifold_tp = mafft_xinsi_plus_centroidalifold_tn = mafft_xinsi_plus_centroidalifold_fp = mafft_xinsi_plus_centroidalifold_fn = 0.
    # ref_sa_plus_centroidalifold_tp = ref_sa_plus_centroidalifold_tn = ref_sa_plus_centroidalifold_fp = ref_sa_plus_centroidalifold_fn = 0.
    for rna_fam_file in os.listdir(rna_fam_dir_path):
      if not rna_fam_file.endswith(".fa"):
        continue
      rna_seq_file_path = os.path.join(rna_fam_dir_path, rna_fam_file)
      rna_seq_lens = [len(rna_seq.seq) for rna_seq in SeqIO.parse(rna_seq_file_path, "fasta")]
      num_of_rnas = len(rna_seq_lens)
      (rna_fam_name, extension) = os.path.splitext(rna_fam_file)
      ref_css_file_path = os.path.join(rna_fam_dir_path, rna_fam_name + ".sth")
      ref_css_and_flat_css = utils.get_css_and_flat_css(utils.get_css_string(ref_css_file_path))
      mearcof_estimated_css_dir_path = os.path.join(mearcof_css_dir_path, "css_of_" + rna_fam_name)
      if not os.path.isdir(mearcof_estimated_css_dir_path):
        continue
      mearcof_estimated_css_dir_path = os.path.join(mearcof_css_dir_path, "css_of_" + rna_fam_name)
      # if not os.path.isdir(centroidalifold_estimated_css_dir_path):
      #   continue
      mearcof_estimated_css_file_path = os.path.join(mearcof_estimated_css_dir_path, "gamma=" + gamma_str + ".dat")
      estimated_css_and_flat_css = utils.get_csss_and_flat_csss(utils.get_css_strings(mearcof_estimated_css_file_path))
      estimated_css, estimated_flat_css = estimated_css_and_flat_css
      ref_css, ref_flat_css = ref_css_and_flat_css
      for m in range(0, num_of_rnas):
        rna_seq_len_1 = rna_seq_lens[m]
        for i in range(0, rna_seq_len_1):
          estimated_bin = (m, i) in estimated_flat_css
          ref_bin = (m, i) in ref_flat_css
          if estimated_bin == ref_bin:
            if estimated_bin == False:
              mearcof_tn += 1
          else:
            if estimated_bin == True:
              mearcof_fp += 1
            else:
              mearcof_fn += 1
          for j in range(i + 1, rna_seq_len_1):
            for n in range(m + 1, num_of_rnas):
              rna_seq_len_2 = rna_seq_lens[n]
              estimated_bin = (m, n, i, j, k, l) in estimated_css
              ref_bin = (m, n, i, j, k, l) in ref_css
              if estimated_bin == ref_bin:
                if estimated_bin == True:
                  mearcof_tp += 1
    ppv = mearcof_tp / (mearcof_tp + mearcof_fp)
    sens = mearcof_tp / (mearcof_tp + mearcof_fn)
    fpr = mearcof_fp / (mearcof_tn + mearcof_fp)
    mearcof_ppvs.insert(0, ppv)
    mearcof_senss.insert(0, sens)
    mearcof_fprs.insert(0, fpr)
    ppv = centroidalifold_tp / (centroidalifold_tp + centroidalifold_fp)
    sens = centroidalifold_tp / (centroidalifold_tp + centroidalifold_fn)
    fpr = centroidalifold_fp / (centroidalifold_tn + centroidalifold_fp)
    centroidalifold_ppvs.insert(0, ppv)
    centroidalifold_senss.insert(0, sens)
    centroidalifold_fprs.insert(0, fpr)
  mearcof_ppvs = numpy.array(mearcof_ppvs) 
  mearcof_senss = numpy.array(mearcof_senss)
  mearcof_fprs = numpy.array(mearcof_fprs)
  # centroidalifold_ppvs = numpy.array(centroidalifold_ppvs) 
  # centroidalifold_senss = numpy.array(centroidalifold_senss)
  # centroidalifold_fprs = numpy.array(centroidalifold_fprs)
  line_1, = pyplot.plot(mearcof_ppvs, mearcof_senss, label = "MEARCOF", marker = "o", linestyle = "-")
  # line_2, = pyplot.plot(centroidalifold_ppvs, centroidalifold_senss, label = "CentroidAlifold", marker = "v", linestyle = "-")
  pyplot.xlabel("Positive predictive value")
  pyplot.ylabel("Sensitivity")
  pyplot.legend(handles = [line_1], loc = 1)
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/ppvs_vs_senss_on_css_estimation.eps", bbox_inches = "tight")
  pyplot.figure()
  line_1, = pyplot.plot(mearcof_fprs, mearcof_senss, label = "MEARCOF", marker = "o", linestyle = "-")
  # line_2, = pyplot.plot(centroidalifold_fprs, centroidalifold_senss, label = "CentroidAlifold", marker = "v", linestyle = "-")
  pyplot.xlabel("False positive rate")
  pyplot.ylabel("Sensitivity")
  pyplot.legend(handles = [line_1], loc = 4)
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/fprs_vs_senss_on_css_estimation.eps", bbox_inches = "tight")

if __name__ == "__main__":
  main()
