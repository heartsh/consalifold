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
  mafft_ginsi_plus_phyloalifold_css_dir_path = asset_dir_path + "/mafft_ginsi_plus_phyloalifold"
  mafft_xinsi_plus_phyloalifold_css_dir_path = asset_dir_path + "/mafft_xinsi_plus_phyloalifold"
  ref_sa_plus_phyloalifold_css_dir_path = asset_dir_path + "/ref_sa_plus_phyloalifold"
  mafft_ginsi_plus_centroidalifold_css_dir_path = asset_dir_path + "/mafft_ginsi_plus_centroidalifold"
  mafft_xinsi_plus_centroidalifold_css_dir_path = asset_dir_path + "/mafft_xinsi_plus_centroidalifold"
  ref_sa_plus_centroidalifold_css_dir_path = asset_dir_path + "/ref_sa_plus_centroidalifold"
  rna_fam_dir_path = asset_dir_path + "/sampled_rna_fams"
  ref_sa_dir_path = asset_dir_path + "/ref_sas"
  mafft_ginsi_plus_phyloalifold_ppvs = []
  mafft_ginsi_plus_phyloalifold_senss = []
  mafft_ginsi_plus_phyloalifold_fprs = []
  mafft_xinsi_plus_phyloalifold_ppvs = []
  mafft_xinsi_plus_phyloalifold_senss = []
  mafft_xinsi_plus_phyloalifold_fprs = []
  ref_sa_plus_phyloalifold_ppvs = []
  ref_sa_plus_phyloalifold_senss = []
  ref_sa_plus_phyloalifold_fprs = []
  mafft_ginsi_plus_centroidalifold_ppvs = []
  mafft_ginsi_plus_centroidalifold_senss = []
  mafft_ginsi_plus_centroidalifold_fprs = []
  mafft_xinsi_plus_centroidalifold_ppvs = []
  mafft_xinsi_plus_centroidalifold_senss = []
  mafft_xinsi_plus_centroidalifold_fprs = []
  ref_sa_plus_centroidalifold_ppvs = []
  ref_sa_plus_centroidalifold_senss = []
  ref_sa_plus_centroidalifold_fprs = []
  gammas = [2. ** i for i in range(-7, 11)]
  for gamma in gammas:
    gamma_str = str(gamma)
    mafft_ginsi_plus_phyloalifold_tp = mafft_ginsi_plus_phyloalifold_tn = mafft_ginsi_plus_phyloalifold_fp = mafft_ginsi_plus_phyloalifold_fn = 0.
    mafft_xinsi_plus_phyloalifold_tp = mafft_xinsi_plus_phyloalifold_tn = mafft_xinsi_plus_phyloalifold_fp = mafft_xinsi_plus_phyloalifold_fn = 0.
    ref_sa_plus_phyloalifold_tp = ref_sa_plus_phyloalifold_tn = ref_sa_plus_phyloalifold_fp = ref_sa_plus_phyloalifold_fn = 0.
    mafft_ginsi_plus_centroidalifold_tp = mafft_ginsi_plus_centroidalifold_tn = mafft_ginsi_plus_centroidalifold_fp = mafft_ginsi_plus_centroidalifold_fn = 0.
    mafft_xinsi_plus_centroidalifold_tp = mafft_xinsi_plus_centroidalifold_tn = mafft_xinsi_plus_centroidalifold_fp = mafft_xinsi_plus_centroidalifold_fn = 0.
    ref_sa_plus_centroidalifold_tp = ref_sa_plus_centroidalifold_tn = ref_sa_plus_centroidalifold_fp = ref_sa_plus_centroidalifold_fn = 0.
    for rna_fam_file in os.listdir(rna_fam_dir_path):
      if not rna_fam_file.endswith(".fa"):
        continue
      rna_seq_file_path = os.path.join(rna_fam_dir_path, rna_fam_file)
      rna_seq_lens = [len(rna_seq.seq) for rna_seq in SeqIO.parse(rna_seq_file_path, "fasta")]
      num_of_rnas = len(rna_seq_lens)
      (rna_fam_name, extension) = os.path.splitext(rna_fam_file)
      ref_css_file_path = os.path.join(ref_sa_dir_path, rna_fam_name + ".sth")
      ref_css_and_flat_css = utils.get_css_and_flat_css(ref_css_file_path)
      mafft_ginsi_plus_phyloalifold_estimated_css_dir_path = os.path.join(mafft_ginsi_plus_phyloalifold_css_dir_path, "csss_of_" + rna_fam_name)
      if not os.path.isdir(mafft_ginsi_plus_phyloalifold_estimated_css_dir_path):
        continue
      mafft_xinsi_plus_phyloalifold_estimated_css_dir_path = os.path.join(mafft_xinsi_plus_phyloalifold_css_dir_path, "csss_of_" + rna_fam_name)
      if not os.path.isdir(mafft_xinsi_plus_phyloalifold_estimated_css_dir_path):
        continue
      ref_sa_plus_phyloalifold_estimated_css_dir_path = os.path.join(ref_sa_plus_phyloalifold_css_dir_path, "csss_of_" + rna_fam_name)
      if not os.path.isdir(ref_sa_plus_phyloalifold_estimated_css_dir_path):
        continue
      mafft_ginsi_plus_centroidalifold_estimated_css_dir_path = os.path.join(mafft_ginsi_plus_centroidalifold_css_dir_path, "csss_of_" + rna_fam_name)
      if not os.path.isdir(mafft_ginsi_plus_centroidalifold_estimated_css_dir_path):
        continue
      mafft_xinsi_plus_centroidalifold_estimated_css_dir_path = os.path.join(mafft_xinsi_plus_centroidalifold_css_dir_path, "csss_of_" + rna_fam_name)
      if not os.path.isdir(mafft_xinsi_plus_centroidalifold_estimated_css_dir_path):
        continue
      ref_sa_plus_centroidalifold_estimated_css_dir_path = os.path.join(ref_sa_plus_centroidalifold_css_dir_path, "csss_of_" + rna_fam_name)
      if not os.path.isdir(ref_sa_plus_centroidalifold_estimated_css_dir_path):
        continue
      mafft_ginsi_plus_phyloalifold_estimated_css_file_path = os.path.join(mafft_ginsi_plus_phyloalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css_and_flat_css = utils.get_css_and_flat_css(mafft_ginsi_plus_phyloalifold_estimated_css_file_path)
      tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
      mafft_ginsi_plus_phyloalifold_tp += tp
      mafft_ginsi_plus_phyloalifold_tn += tn
      mafft_ginsi_plus_phyloalifold_fp += fp
      mafft_ginsi_plus_phyloalifold_fn += fn
      mafft_xinsi_plus_phyloalifold_estimated_css_file_path = os.path.join(mafft_xinsi_plus_phyloalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css_and_flat_css = utils.get_css_and_flat_css(mafft_xinsi_plus_phyloalifold_estimated_css_file_path)
      tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
      mafft_xinsi_plus_phyloalifold_tp += tp
      mafft_xinsi_plus_phyloalifold_tn += tn
      mafft_xinsi_plus_phyloalifold_fp += fp
      mafft_xinsi_plus_phyloalifold_fn += fn
      ref_sa_plus_phyloalifold_estimated_css_file_path = os.path.join(ref_sa_plus_phyloalifold_estimated_css_dir_path, "gamma=" + gamma_str + ".sth")
      estimated_css_and_flat_css = utils.get_css_and_flat_css(ref_sa_plus_phyloalifold_estimated_css_file_path)
      tp, tn, fp, fn = get_bin_counts(rna_seq_lens, estimated_css_and_flat_css, ref_css_and_flat_css)
      ref_sa_plus_phyloalifold_tp += tp
      ref_sa_plus_phyloalifold_tn += tn
      ref_sa_plus_phyloalifold_fp += fp
      ref_sa_plus_phyloalifold_fn += fn
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
    ppv = mafft_ginsi_plus_phyloalifold_tp / (mafft_ginsi_plus_phyloalifold_tp + mafft_ginsi_plus_phyloalifold_fp)
    sens = mafft_ginsi_plus_phyloalifold_tp / (mafft_ginsi_plus_phyloalifold_tp + mafft_ginsi_plus_phyloalifold_fn)
    fpr = mafft_ginsi_plus_phyloalifold_fp / (mafft_ginsi_plus_phyloalifold_tn + mafft_ginsi_plus_phyloalifold_fp)
    mafft_ginsi_plus_phyloalifold_ppvs.insert(0, ppv)
    mafft_ginsi_plus_phyloalifold_senss.insert(0, sens)
    mafft_ginsi_plus_phyloalifold_fprs.insert(0, fpr)
    ppv = mafft_xinsi_plus_phyloalifold_tp / (mafft_xinsi_plus_phyloalifold_tp + mafft_xinsi_plus_phyloalifold_fp)
    sens = mafft_xinsi_plus_phyloalifold_tp / (mafft_xinsi_plus_phyloalifold_tp + mafft_xinsi_plus_phyloalifold_fn)
    fpr = mafft_xinsi_plus_phyloalifold_fp / (mafft_xinsi_plus_phyloalifold_tn + mafft_xinsi_plus_phyloalifold_fp)
    mafft_xinsi_plus_phyloalifold_ppvs.insert(0, ppv)
    mafft_xinsi_plus_phyloalifold_senss.insert(0, sens)
    mafft_xinsi_plus_phyloalifold_fprs.insert(0, fpr)
    ppv = ref_sa_plus_phyloalifold_tp / (ref_sa_plus_phyloalifold_tp + ref_sa_plus_phyloalifold_fp)
    sens = ref_sa_plus_phyloalifold_tp / (ref_sa_plus_phyloalifold_tp + ref_sa_plus_phyloalifold_fn)
    fpr = ref_sa_plus_phyloalifold_fp / (ref_sa_plus_phyloalifold_tn + ref_sa_plus_phyloalifold_fp)
    ref_sa_plus_phyloalifold_ppvs.insert(0, ppv)
    ref_sa_plus_phyloalifold_senss.insert(0, sens)
    ref_sa_plus_phyloalifold_fprs.insert(0, fpr)
    ppv = mafft_ginsi_plus_centroidalifold_tp / (mafft_ginsi_plus_centroidalifold_tp + mafft_ginsi_plus_centroidalifold_fp)
    sens = mafft_ginsi_plus_centroidalifold_tp / (mafft_ginsi_plus_centroidalifold_tp + mafft_ginsi_plus_centroidalifold_fn)
    fpr = mafft_ginsi_plus_centroidalifold_fp / (mafft_ginsi_plus_centroidalifold_tn + mafft_ginsi_plus_centroidalifold_fp)
    mafft_ginsi_plus_centroidalifold_ppvs.insert(0, ppv)
    mafft_ginsi_plus_centroidalifold_senss.insert(0, sens)
    mafft_ginsi_plus_centroidalifold_fprs.insert(0, fpr)
    ppv = mafft_xinsi_plus_centroidalifold_tp / (mafft_xinsi_plus_centroidalifold_tp + mafft_xinsi_plus_centroidalifold_fp)
    sens = mafft_xinsi_plus_centroidalifold_tp / (mafft_xinsi_plus_centroidalifold_tp + mafft_xinsi_plus_centroidalifold_fn)
    fpr = mafft_xinsi_plus_centroidalifold_fp / (mafft_xinsi_plus_centroidalifold_tn + mafft_xinsi_plus_centroidalifold_fp)
    mafft_xinsi_plus_centroidalifold_ppvs.insert(0, ppv)
    mafft_xinsi_plus_centroidalifold_senss.insert(0, sens)
    mafft_xinsi_plus_centroidalifold_fprs.insert(0, fpr)
    ppv = ref_sa_plus_centroidalifold_tp / (ref_sa_plus_centroidalifold_tp + ref_sa_plus_centroidalifold_fp)
    sens = ref_sa_plus_centroidalifold_tp / (ref_sa_plus_centroidalifold_tp + ref_sa_plus_centroidalifold_fn)
    fpr = ref_sa_plus_centroidalifold_fp / (ref_sa_plus_centroidalifold_tn + ref_sa_plus_centroidalifold_fp)
    ref_sa_plus_centroidalifold_ppvs.insert(0, ppv)
    ref_sa_plus_centroidalifold_senss.insert(0, sens)
    ref_sa_plus_centroidalifold_fprs.insert(0, fpr)
  mafft_ginsi_plus_phyloalifold_ppvs = numpy.array(mafft_ginsi_plus_phyloalifold_ppvs)
  mafft_ginsi_plus_phyloalifold_senss = numpy.array(mafft_ginsi_plus_phyloalifold_senss)
  mafft_ginsi_plus_phyloalifold_fprs = numpy.array(mafft_ginsi_plus_phyloalifold_fprs)
  mafft_xinsi_plus_phyloalifold_ppvs = numpy.array(mafft_xinsi_plus_phyloalifold_ppvs)
  mafft_xinsi_plus_phyloalifold_senss = numpy.array(mafft_xinsi_plus_phyloalifold_senss)
  mafft_xinsi_plus_phyloalifold_fprs = numpy.array(mafft_xinsi_plus_phyloalifold_fprs)
  ref_sa_plus_phyloalifold_ppvs = numpy.array(ref_sa_plus_phyloalifold_ppvs)
  ref_sa_plus_phyloalifold_senss = numpy.array(ref_sa_plus_phyloalifold_senss)
  ref_sa_plus_phyloalifold_fprs = numpy.array(ref_sa_plus_phyloalifold_fprs)
  mafft_ginsi_plus_centroidalifold_ppvs = numpy.array(mafft_ginsi_plus_centroidalifold_ppvs)
  mafft_ginsi_plus_centroidalifold_senss = numpy.array(mafft_ginsi_plus_centroidalifold_senss)
  mafft_ginsi_plus_centroidalifold_fprs = numpy.array(mafft_ginsi_plus_centroidalifold_fprs)
  mafft_xinsi_plus_centroidalifold_ppvs = numpy.array(mafft_xinsi_plus_centroidalifold_ppvs)
  mafft_xinsi_plus_centroidalifold_senss = numpy.array(mafft_xinsi_plus_centroidalifold_senss)
  mafft_xinsi_plus_centroidalifold_fprs = numpy.array(mafft_xinsi_plus_centroidalifold_fprs)
  ref_sa_plus_centroidalifold_ppvs = numpy.array(ref_sa_plus_centroidalifold_ppvs)
  ref_sa_plus_centroidalifold_senss = numpy.array(ref_sa_plus_centroidalifold_senss)
  ref_sa_plus_centroidalifold_fprs = numpy.array(ref_sa_plus_centroidalifold_fprs)
  line_1, = pyplot.plot(ref_sa_plus_phyloalifold_ppvs, ref_sa_plus_phyloalifold_senss, label = "Reference + PhyloAliFold", marker = "^", linestyle = "-")
  line_2, = pyplot.plot(mafft_xinsi_plus_phyloalifold_ppvs, mafft_xinsi_plus_phyloalifold_senss, label = "MAFFT X-INS-i + PhyloAliFold", marker = "v", linestyle = "-")
  line_3, = pyplot.plot(mafft_ginsi_plus_phyloalifold_ppvs, mafft_ginsi_plus_phyloalifold_senss, label = "MAFFT G-INS-i + PhyloAliFold", marker = "o", linestyle = "-")
  line_4, = pyplot.plot(ref_sa_plus_centroidalifold_ppvs, ref_sa_plus_centroidalifold_senss, label = "Reference + CentroidAlifold", marker = "D", linestyle = "-")
  line_5, = pyplot.plot(mafft_xinsi_plus_centroidalifold_ppvs, mafft_xinsi_plus_centroidalifold_senss, label = "MAFFT X-INS-i + CentroidAlifold", marker = "p", linestyle = "-")
  line_6, = pyplot.plot(mafft_ginsi_plus_centroidalifold_ppvs, mafft_ginsi_plus_centroidalifold_senss, label = "MAFFT G-INS-i + CentroidAlifold", marker = "s", linestyle = "-")
  pyplot.xlabel("Positive predictive value")
  pyplot.ylabel("Sensitivity")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4, line_5, line_6], loc = 3)
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  pyplot.tight_layout()
  pyplot.savefig(image_dir_path + "/ppvs_vs_senss_on_css_estimation.eps", bbox_inches = "tight")
  pyplot.figure()
  line_1, = pyplot.plot(ref_sa_plus_phyloalifold_fprs, ref_sa_plus_phyloalifold_senss, label = "Reference + PhyloAliFold", marker = "^", linestyle = "-")
  line_2, = pyplot.plot(mafft_xinsi_plus_phyloalifold_fprs, mafft_xinsi_plus_phyloalifold_senss, label = "MAFFT X-INS-i + PhyloAliFold", marker = "v", linestyle = "-")
  line_3, = pyplot.plot(mafft_ginsi_plus_phyloalifold_fprs, mafft_ginsi_plus_phyloalifold_senss, label = "MAFFT G-INS-i + PhyloAliFold", marker = "o", linestyle = "-")
  line_4, = pyplot.plot(ref_sa_plus_centroidalifold_fprs, ref_sa_plus_centroidalifold_senss, label = "Reference + CentroidAlifold", marker = "D", linestyle = "-")
  line_5, = pyplot.plot(mafft_xinsi_plus_centroidalifold_fprs, mafft_xinsi_plus_centroidalifold_senss, label = "MAFFT X-INS-i + CentroidAlifold", marker = "p", linestyle = "-")
  line_6, = pyplot.plot(mafft_ginsi_plus_centroidalifold_fprs, mafft_ginsi_plus_centroidalifold_senss, label = "MAFFT G-INS-i + CentroidAlifold", marker = "s", linestyle = "-")
  pyplot.xlabel("False positive rate")
  pyplot.ylabel("Sensitivity")
  pyplot.legend(handles = [line_1, line_2, line_3, line_4, line_5, line_6], loc = 4)
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
