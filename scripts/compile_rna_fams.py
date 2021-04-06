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
from Bio.Align import MultipleSeqAlignment
from sklearn.model_selection import train_test_split

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  rfam_seed_sta_file_path = asset_dir_path + "/rfam_seed_stas_v14.3.sth"
  rna_seq_dir_path_valid = asset_dir_path + "/compiled_rna_fams_valid"
  rna_seq_dir_path_test = asset_dir_path + "/compiled_rna_fams_test"
  ref_sa_dir_path_valid = asset_dir_path + "/ref_sas_valid"
  ref_sa_dir_path_test = asset_dir_path + "/ref_sas_test"
  if not os.path.isdir(rna_seq_dir_path_valid):
    os.mkdir(rna_seq_dir_path_valid)
  if not os.path.isdir(rna_seq_dir_path_test):
    os.mkdir(rna_seq_dir_path_test)
  if not os.path.isdir(ref_sa_dir_path_valid):
    os.mkdir(ref_sa_dir_path_valid)
  if not os.path.isdir(ref_sa_dir_path_test):
    os.mkdir(ref_sa_dir_path_test)
  max_sa_len = 500
  max_seq_num = 20
  stas = [sta for sta in AlignIO.parse(rfam_seed_sta_file_path, "stockholm") if len(sta[0]) <= max_sa_len and len(sta) <= max_seq_num and is_valid(sta)]
  num_of_stas = len(stas)
  print("# RNA families: %d" % num_of_stas)
  valid_data, test_data = train_test_split(stas, test_size = 0.5)
  for i, sta in enumerate(valid_data):
    sa_file_path = os.path.join(ref_sa_dir_path_valid, "rna_fam_%d.sth" % i)
    AlignIO.write(sta, sa_file_path, "stockholm")
    sa_file_path = os.path.join(ref_sa_dir_path_valid, "rna_fam_%d.fa" % i)
    AlignIO.write(sta, sa_file_path, "fasta")
    sa_file_path = os.path.join(ref_sa_dir_path_valid, "rna_fam_%d.aln" % i)
    AlignIO.write(sta, sa_file_path, "clustal")
    rna_seq_file_path = os.path.join(rna_seq_dir_path_valid, "rna_fam_%d.fa" % i)
    rna_seq_file = open(rna_seq_file_path, "w")
    for j, rec in enumerate(sta):
      rna_seq_file.write(">%d(%s)\n%s\n" % (j, rec.id, str(rec.seq).replace("-", "")))
  for i, sta in enumerate(test_data):
    sa_file_path = os.path.join(ref_sa_dir_path_test, "rna_fam_%d.sth" % i)
    AlignIO.write(sta, sa_file_path, "stockholm")
    sa_file_path = os.path.join(ref_sa_dir_path_test, "rna_fam_%d.fa" % i)
    AlignIO.write(sta, sa_file_path, "fasta")
    sa_file_path = os.path.join(ref_sa_dir_path_test, "rna_fam_%d.aln" % i)
    AlignIO.write(sta, sa_file_path, "clustal")
    rna_seq_file_path = os.path.join(rna_seq_dir_path_test, "rna_fam_%d.fa" % i)
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
