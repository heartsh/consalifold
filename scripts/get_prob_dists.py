#! /usr/bin/env python

import utils
from Bio import SeqIO
import numpy
import seaborn
from matplotlib import pyplot
import matplotlib
import os
import multiprocessing
import time
import datetime
import shutil

node_size = 150
edge_width = 2
label_font_size = 10
edge_label_font_size = 4
cmap = pyplot.cm.viridis
seaborn.set()
color_palette = seaborn.color_palette()
white = "#F2F2F2"

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  seq_file_path = asset_dir_path + "/test_seqs.fa"
  seqs = [rec for rec in SeqIO.parse(seq_file_path, "fasta")]
  seq = seqs[0]
  seq_lens = [len(seq) for seq in seqs]
  seq_len = len(seq)
  upp_mat_4_hl_turner = utils.get_upp_mats(asset_dir_path + "/test_seqs/upp_mats_on_hl.dat", seq_lens)[0]
  bpp_mat_2_turner = utils.get_upp_mats(asset_dir_path + "/test_seqs/bpp_mats_2.dat", seq_lens)[0]
  upp_mat_4_bl_turner = utils.get_upp_mats(asset_dir_path + "/test_seqs/upp_mats_on_bl.dat", seq_lens)[0]
  upp_mat_4_il_turner = utils.get_upp_mats(asset_dir_path + "/test_seqs/upp_mats_on_il.dat", seq_lens)[0]
  upp_mat_4_ml_turner = utils.get_upp_mats(asset_dir_path + "/test_seqs/upp_mats_on_ml.dat", seq_lens)[0]
  upp_mat_4_el_turner = utils.get_upp_mats(asset_dir_path + "/test_seqs/upp_mats_on_el.dat", seq_lens)[0]
  pyplot.figure(figsize=(7, 7))
  pyplot.stackplot(range(seq_len), upp_mat_4_hl_turner, bpp_mat_2_turner, upp_mat_4_bl_turner, upp_mat_4_il_turner, upp_mat_4_ml_turner, upp_mat_4_el_turner)
  legends = ["Unpairing in hairpin loop", "Base-pairing", "Unpairing in bulge loop", "Unpairing in interior loop", "Unpairing in multi-loop", "Unpairing in external loop"]
  pyplot.legend(legends, loc = "upper right", bbox_to_anchor=(1.0, 1.23))
  pyplot.savefig(image_dir_path + "/consprob_struct_context_profs_turner.png", bbox_inches = "tight")
  pyplot.clf()
  capr_prof_seqs = utils.get_capr_prof_seqs(asset_dir_path + "/capr_test_seqs.dat")
  pyplot.figure(figsize=(11.5, 11.5))
  pyplot.stackplot(range(seq_len), capr_prof_seqs["Hairpin"], capr_prof_seqs["Stem"], capr_prof_seqs["Bulge"], capr_prof_seqs["Internal"], capr_prof_seqs["Multibranch"], capr_prof_seqs["Exterior"])
  legends = ["Unpairing in hairpin loop", "Base-pairing", "Unpairing in bulge loop", "Unpairing in interior loop", "Unpairing in multi-loop", "Unpairing in external loop"]
  pyplot.legend(legends, loc = "upper right", bbox_to_anchor=(1.0, 1.13))
  pyplot.savefig(image_dir_path + "/capr_struct_context_profs.png", bbox_inches = "tight")

if __name__ == "__main__":
  main()
