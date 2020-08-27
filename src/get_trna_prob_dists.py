#! /usr/bin/env python

import utils
from Bio import SeqIO
import seaborn
from matplotlib import pyplot
import os

seaborn.set()
cmap = pyplot.cm.viridis

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  seq_file_path = asset_dir_path + "/sampled_trnas.fa"
  seqs = [rec for rec in SeqIO.parse(seq_file_path, "fasta")]
  seq = seqs[0]
  seq_lens = [len(seq) for seq in seqs]
  seq_len = len(seq)
  bpp_mat = utils.get_bpp_mats(asset_dir_path + "/sampled_trnas/bpp_mats.dat", seq_lens)[0]
  bpp_mat_on_ss = utils.get_bpp_mats(asset_dir_path + "/sampled_trnas/bpp_mats_on_ss.dat", seq_lens)[0]
  xlabels = [str(i + 1) if (i + 1) % 20 == 0 else "" for i in range(seq_len)]
  (_, axes) = pyplot.subplots(nrows = 1, ncols = 2, figsize = (12, 6))
  seaborn.heatmap(bpp_mat_on_ss, ax = axes[0], xticklabels = xlabels, yticklabels = xlabels, cbar = False, cmap = cmap, vmin = 0, vmax = 1)
  seaborn.heatmap(bpp_mat, ax = axes[1], xticklabels = xlabels, yticklabels = xlabels, cbar = False, cmap = cmap, vmin = 0, vmax = 1)
  pyplot.savefig(image_dir_path + "/trna_bpp_mat_comparison.eps", bbox_inches = "tight")
  pyplot.clf()

if __name__ == "__main__":
  main()
