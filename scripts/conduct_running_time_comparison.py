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
from scipy import stats
import pandas

seaborn.set()

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  num_of_threads = multiprocessing.cpu_count()
  gammas = [2. ** i for i in range(-7, 11)]
  mafft_xinsi_params = []
  consalifold_params = []
  posterior_consalifold_params = []
  rna_seq_dir_path = asset_dir_path + "/compiled_rna_fams_test"
  mafft_xinsi_dir_path = asset_dir_path + "/mafft_xinsi"
  mafft_xinsi_plus_consalifold_dir_path = asset_dir_path + "/mafft_xinsi_plus_consalifold"
  posterior_mafft_xinsi_plus_consalifold_dir_path = asset_dir_path + "/posterior_mafft_xinsi_plus_consalifold"
  if not os.path.isdir(mafft_xinsi_plus_consalifold_dir_path):
    os.mkdir(mafft_xinsi_plus_consalifold_dir_path)
  if not os.path.isdir(posterior_mafft_xinsi_plus_consalifold_dir_path):
    os.mkdir(posterior_mafft_xinsi_plus_consalifold_dir_path)
  sub_thread_num = 4
  for rna_seq_file in os.listdir(rna_seq_dir_path):
    if not rna_seq_file.endswith(".fa"):
      continue
    rna_seq_file_path = os.path.join(rna_seq_dir_path, rna_seq_file)
    (rna_family_name, extension) = os.path.splitext(rna_seq_file)
    mafft_xinsi_output_file_path = os.path.join(mafft_xinsi_dir_path, rna_family_name + ".aln")
    mafft_xinsi_plus_consalifold_output_dir_path = os.path.join(mafft_xinsi_plus_consalifold_dir_path, rna_family_name)
    posterior_mafft_xinsi_plus_consalifold_output_dir_path = os.path.join(posterior_mafft_xinsi_plus_consalifold_dir_path, rna_family_name)
    if not os.path.isdir(mafft_xinsi_plus_consalifold_output_dir_path):
      os.mkdir(mafft_xinsi_plus_consalifold_output_dir_path)
    if not os.path.isdir(posterior_mafft_xinsi_plus_consalifold_output_dir_path):
      os.mkdir(posterior_mafft_xinsi_plus_consalifold_output_dir_path)
    consalifold_params.insert(0, (sub_thread_num, mafft_xinsi_output_file_path, mafft_xinsi_plus_consalifold_output_dir_path, False))
    posterior_consalifold_params.insert(0, (sub_thread_num, mafft_xinsi_output_file_path, posterior_mafft_xinsi_plus_consalifold_output_dir_path, True))
  # ConsAliFold's execution.
  pool = multiprocessing.Pool(int(num_of_threads / sub_thread_num))
  consalifold_results = pool.map(bench_consalifold, consalifold_params)
  consalifold_output_file_path = asset_dir_path + "/consalifold_running_times_turner.dat"
  write_consalifold_results(consalifold_results, consalifold_output_file_path)
  data_turner = read_consalifold_results(consalifold_output_file_path)
  posterior_consalifold_results = pool.map(bench_consalifold, posterior_consalifold_params)
  posterior_consalifold_output_file_path = asset_dir_path + "/consalifold_running_times_posterior.dat"
  write_consalifold_results(posterior_consalifold_results, posterior_consalifold_output_file_path)
  data_posterior = read_consalifold_results(posterior_consalifold_output_file_path)
  data = {"Running time (s)": data_turner + data_posterior, "Pair-matching probability inference method": ["ConsProb"] * len(data_turner) + ["LocARNA-P + our PCT"] * len(data_posterior)}
  data_frame = pandas.DataFrame(data = data)
  ax = seaborn.boxplot(x = "Pair-matching probability inference method", y = "Running time (s)", data = data_frame, sym = "")
  fig = ax.get_figure()
  fig.tight_layout()
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  fig.savefig(image_dir_path + "/consalifold_model_comparison_running_time.eps", bbox_inches = "tight")
  fig.clf()
  print("Running time significance test: ", stats.ttest_rel(data_turner, data_posterior))

def read_consalifold_results(consalifold_output_file_path):
  with open(consalifold_output_file_path) as f:
    lines = f.readlines()
    lines = [float(line) for line in lines]
    return lines

def write_consalifold_results(consalifold_results, consalifold_output_file_path):
  with open(consalifold_output_file_path, "w") as f:
    buf = ""
    for consalifold_result in consalifold_results:
      buf += "%f\n" % consalifold_result
    f.write(buf)

def bench_consalifold(consalifold_params):
  (sub_thread_num, input_sa_file_path, consalifold_output_dir_path, is_posterior_model) = consalifold_params
  consalifold_command = "consalifold %s-b -t " % ("-u " if is_posterior_model else "") + str(sub_thread_num) + " -i " + input_sa_file_path + " -o " + consalifold_output_dir_path
  begin = time.time()
  utils.run_command(consalifold_command)
  consalifold_elapsed_time = time.time() - begin
  return consalifold_elapsed_time

def run_mafft_xinsi(mafft_xinsi_params):
  (rna_seq_file_path, mafft_xinsi_output_file_path) = mafft_xinsi_params
  mafft_xinsi_command = "mafft-xinsi --thread 1 --quiet " + "--clustalout " + rna_seq_file_path + " > " + mafft_xinsi_output_file_path
  utils.run_command(mafft_xinsi_command)
  sa = AlignIO.read(mafft_xinsi_output_file_path, "clustal")
  mafft_xinsi_output_file_path = os.path.splitext(mafft_xinsi_output_file_path)[0] + ".fa"
  AlignIO.write(sa, mafft_xinsi_output_file_path, "fasta")

if __name__ == "__main__":
  main()
