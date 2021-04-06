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

taus = numpy.arange(0, 1.1, 0.1).tolist()

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  num_of_threads = multiprocessing.cpu_count()
  consalifold_params = []
  rna_seq_dir_path = asset_dir_path + "/compiled_rna_fams"
  # rna_seq_dir_path = asset_dir_path + "/compiled_rna_fams_4_micro_bench"
  ref_sa_dir_path = asset_dir_path + "/ref_sas"
  consalifold_dir_path = asset_dir_path + "/grid_search"
  if not os.path.isdir(consalifold_dir_path):
    os.mkdir(consalifold_dir_path)
  sub_thread_num = 4
  for tau in taus:
    tau_str = str(tau)
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

if __name__ == "__main__":
  main()
