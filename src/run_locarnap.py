#! /usr/bin/env python

import utils
from Bio import SeqIO
import seaborn
from matplotlib import pyplot
import os
import multiprocessing
import time
import datetime
import shutil
from os import path

seaborn.set()
pyplot.rcParams['legend.handlelength'] = 0

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  num_of_threads = multiprocessing.cpu_count()
  locarnap_output_dir_path = asset_dir_path + "/locarnap"
  if not os.path.isdir(locarnap_output_dir_path):
    os.mkdir(locarnap_output_dir_path)
  temp_dir_path = "/tmp/run_locarnap_%s" % datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')
  if not os.path.isdir(temp_dir_path):
    os.mkdir(temp_dir_path)
  rna_seq_dir_path = asset_dir_path + "/compiled_rna_fams"
  locarnap_params = []
  for rna_seq_file in os.listdir(rna_seq_dir_path):
    if not rna_seq_file.endswith(".fa"):
      continue
    rna_seq_file_path = os.path.join(rna_seq_dir_path, rna_seq_file)
    (rna_family_name, extension) = os.path.splitext(rna_seq_file)
    locarnap_output_file_path = os.path.join(locarnap_output_dir_path, rna_family_name + "_bpp_mats.dat")
    locarnap_params.insert(0, (rna_seq_file_path, locarnap_output_file_path, temp_dir_path))
  pool = multiprocessing.Pool(num_of_threads)
  pool.map(run_locarnap, locarnap_params)
  shutil.rmtree(temp_dir_path)

def run_locarnap(locarnap_params):
  (seq_file_path, locarnap_output_file_path, temp_dir_path) = locarnap_params
  base_name = os.path.basename(seq_file_path)
  rna_family = os.path.splitext(base_name)[0]
  sub_temp_dir_path = temp_dir_path + "/" + rna_family
  if not os.path.isdir(sub_temp_dir_path):
    os.mkdir(sub_temp_dir_path)
  recs = [rec for rec in SeqIO.parse(seq_file_path, "fasta")]
  num_of_recs = len(recs)
  bpp_mats = []
  for i in range(num_of_recs):
    rec_1 = recs[i]
    seq_file_path = os.path.join(sub_temp_dir_path, "seq_%d.fa" % i)
    SeqIO.write(rec_1, open(seq_file_path, "w"), "fasta")
    bpp_mats.append({})
  for i in range(num_of_recs):
    seq_file_path_1 = os.path.join(sub_temp_dir_path, "seq_%d.fa" % i)
    bpp_mat = bpp_mats[i]
    for j in range(i + 1, num_of_recs):
      bpp_mat_2 = bpp_mats[j]
      seq_file_path_2 = os.path.join(sub_temp_dir_path, "seq_%d.fa" % j)
      temp_locarnap_output_file_path = os.path.join(sub_temp_dir_path, "locarnap_output_%d%d.dat" % (i, j))
      locarnap_command = "locarna_p " + seq_file_path_1 + " " + seq_file_path_2 + " --write-arcmatch-probs=" + temp_locarnap_output_file_path
      utils.run_command(locarnap_command)
      bpap_mat = read_locarnap_output(temp_locarnap_output_file_path)
      for ((m, n, o, p), bpap) in bpap_mat.items():
        term = bpap / (num_of_recs - 1)
        if True:
          if (m, n) in bpp_mat:
            bpp_mat[(m, n)] += term
          else:
            bpp_mat[(m, n)] = term
          if (o, p) in bpp_mat_2:
            bpp_mat_2[(o, p)] += term
          else:
            bpp_mat_2[(o, p)] = term
      bpp_mats[j] = bpp_mat_2
    bpp_mats[i] = bpp_mat
  locarnap_output_file = open(locarnap_output_file_path, "w")
  contents = ""
  for (i) in range(num_of_recs):
    bpp_mat = bpp_mats[i]
    contents += ">%d\n" % i if i == 0 else "\n\n>%d\n" % i
    for (j, ((m, n), bpp)) in enumerate(bpp_mat.items()):
      contents += "%d,%d,%f" % (m, n, bpp) if j == 0 else " %d,%d,%f" % (m, n, bpp)
  locarnap_output_file.write(contents)
  locarnap_output_file.close()

def read_locarnap_output(locarnap_output_file_path):
  bpap_mat = dict()
  locarnap_output_file = open(locarnap_output_file_path)
  for line in locarnap_output_file.readlines():
    split = line.split()
    i, j, k, l, bpap = int(split[0]), int(split[1]), int(split[2]), int(split[3]), float(split[4])
    bpap_mat[(i, j, k, l)] = bpap
  locarnap_output_file.close()
  return bpap_mat

if __name__ == "__main__":
  main()
