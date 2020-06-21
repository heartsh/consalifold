#! /usr/bin/env python

import utils
from Bio import SeqIO
from Bio import AlignIO
import numpy
import seaborn
from matplotlib import pyplot
import os
import multiprocessing
import time
import datetime
import shutil

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  rfam_seed_sta_file_path = asset_dir_path + "/rfam_seed_stas.sth"
  for data_set in ["short", "long"]:
    rna_seq_dir_path = asset_dir_path + "/compiled_rna_fams_" + data_set
    rna_seq_dir_path_4_micro_bench = asset_dir_path + "/compiled_rna_fams_" + data_set + "_4_micro_bench"
    ref_sa_dir_path = asset_dir_path + "/ref_sas_" + data_set
    ref_sa_dir_path_4_micro_bench = asset_dir_path + "/ref_sas_" + data_set + "_4_micro_bench"
    if not os.path.isdir(rna_seq_dir_path):
      os.mkdir(rna_seq_dir_path)
    if not os.path.isdir(rna_seq_dir_path_4_micro_bench):
      os.mkdir(rna_seq_dir_path_4_micro_bench)
    if not os.path.isdir(ref_sa_dir_path):
      os.mkdir(ref_sa_dir_path)
    if not os.path.isdir(ref_sa_dir_path_4_micro_bench):
      os.mkdir(ref_sa_dir_path_4_micro_bench)
    min_sa_len = 0 if data_set == "short" else 201
    max_sa_len = 200 if data_set == "short" else 500
    max_seq_num = 10 if data_set == "short" else 5
    stas = [sta for sta in AlignIO.parse(rfam_seed_sta_file_path, "stockholm") if len(sta) <= max_seq_num and min_sa_len <= len(sta[0]) and len(sta[0]) <= max_sa_len and is_valid(sta)]
    num_of_stas = len(stas)
    print("# %s RNA families: %d" % (data_set, num_of_stas))
    sample_rate = 0.05 if data_set == "short" else 0.1
    num_of_samples = int(sample_rate * num_of_stas)
    print("# %s RNA families for micro benchmark: %d" % (data_set, num_of_samples))
    sampled_stas = numpy.random.choice(stas, num_of_samples, replace = False)
    for i, sta in enumerate(stas):
      sa_file_path = os.path.join(ref_sa_dir_path, "rna_fam_%d.sth" % i)
      AlignIO.write(sta, sa_file_path, "stockholm")
      sa_file_path = os.path.join(ref_sa_dir_path, "rna_fam_%d.aln" % i)
      AlignIO.write(sta, sa_file_path, "clustal")
      rna_seq_file_path = os.path.join(rna_seq_dir_path, "rna_fam_%d.fa" % i)
      rna_seq_file = open(rna_seq_file_path, "w")
      for j, rec in enumerate(sta):
        rna_seq_file.write(">%d(%s)\n%s\n" % (j, rec.id, str(rec.seq).replace("-", "")))
    for i, sta in enumerate(sampled_stas):
      sa_file_path = os.path.join(ref_sa_dir_path_4_micro_bench, "rna_fam_%d.sth" % i)
      AlignIO.write(sta, sa_file_path, "stockholm")
      sa_file_path = os.path.join(ref_sa_dir_path_4_micro_bench, "rna_fam_%d.aln" % i)
      AlignIO.write(sta, sa_file_path, "clustal")
      rna_seq_file_path = os.path.join(rna_seq_dir_path_4_micro_bench, "rna_fam_%d.fa" % i)
      rna_seq_file = open(rna_seq_file_path, "w")
      for j, rec in enumerate(sta):
        rna_seq_file.write(">%d(%s)\n%s\n" % (j, rec.id, str(rec.seq).replace("-", "")))

def is_valid(sta):
  for row in sta:
    if any(char in str(row.seq) for char in "RYWSMKHBVDN"):
      return False
  return True

if __name__ == "__main__":
  main()
