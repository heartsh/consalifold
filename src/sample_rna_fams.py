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
  rna_seq_dir_path = asset_dir_path + "/sampled_rna_fams"
  ref_sa_dir_path = asset_dir_path + "/ref_sas"
  if not os.path.isdir(rna_seq_dir_path):
    os.mkdir(rna_seq_dir_path)
  if not os.path.isdir(ref_sa_dir_path):
    os.mkdir(ref_sa_dir_path)
  rfam_seed_sta_file_path = asset_dir_path + "/rfam_seed_stas.sth"
  max_seq_num = 10
  max_sa_len = 200
  stas = [sta for sta in AlignIO.parse(rfam_seed_sta_file_path, "stockholm") if len(sta) <= max_seq_num and len(sta[0]) <= max_sa_len and is_valid(sta)]
  num_of_stas = len(stas)
  print(num_of_stas)
  sample_rate = 0.1
  num_of_samples = int(sample_rate * num_of_stas)
  print(num_of_samples)
  sampled_stas = numpy.random.choice(stas, num_of_samples, replace = False)
  for i, sta in enumerate(sampled_stas):
    sa_file_path = os.path.join(ref_sa_dir_path, "rna_fam_%d.sth" % i)
    AlignIO.write(sta, sa_file_path, "stockholm")
    sa_file_path = os.path.join(ref_sa_dir_path, "rna_fam_%d.aln" % i)
    AlignIO.write(sta, sa_file_path, "clustal")
    rna_seq_file_path = os.path.join(rna_seq_dir_path, "rna_fam_%d.fa" % i)
    rna_seq_file = open(rna_seq_file_path, "w")
    for i, rec in enumerate(sta):
      rna_seq_file.write(">%d(%s)\n%s\n" % (i, rec.id, str(rec.seq).replace("-", "")))

def is_valid(sta):
  for row in sta:
    if any(char in str(row.seq) for char in "RYWSMKHBVDN"):
      return False
  return True

if __name__ == "__main__":
  main()
