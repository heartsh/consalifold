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
  consalifold_params = []
  contra_consalifold_params = []
  ref_sa_dir_path = asset_dir_path + "/ref_sas_all"
  ref_sa_plus_consalifold_dir_path = asset_dir_path + "/ref_sa_plus_consalifold_all"
  if not os.path.isdir(ref_sa_plus_consalifold_dir_path):
    os.mkdir(ref_sa_plus_consalifold_dir_path)
  sub_thread_num = 4
  for ref_sa_file in os.listdir(ref_sa_dir_path):
    if not ref_sa_file.endswith(".aln"):
      continue
    ref_sa_file_path = os.path.join(ref_sa_dir_path, ref_sa_file)
    (rna_family_name, extension) = os.path.splitext(ref_sa_file)
    ref_sa_plus_consalifold_output_dir_path = os.path.join(ref_sa_plus_consalifold_dir_path, rna_family_name)
    if not os.path.isdir(ref_sa_plus_consalifold_output_dir_path):
      os.mkdir(ref_sa_plus_consalifold_output_dir_path)
    consalifold_params.insert(0, (sub_thread_num, ref_sa_file_path, ref_sa_plus_consalifold_output_dir_path, False))
  # ConsAliFold's execution.
  pool = multiprocessing.Pool(int(num_of_threads / sub_thread_num))
  consalifold_results = pool.map(bench_consalifold, consalifold_params)
  consalifold_output_file_path = asset_dir_path + "/consalifold_running_times_turner_3.dat"
  write_consalifold_results(consalifold_results, consalifold_output_file_path)
  data = read_consalifold_results(consalifold_output_file_path)
  data = {"Running time (s)": data[0], "Maximum sequence length": data[1], "Number of RNA sequences": data[2]}
  data_frame = pandas.DataFrame(data = data)
  ax = seaborn.scatterplot(data = data_frame, x = "Maximum sequence length", y = "Running time (s)", markers = True)
  fig = ax.get_figure()
  fig.tight_layout()
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  fig.savefig(image_dir_path + "/consalifold_running_time_max_seq_length.eps", bbox_inches = "tight")
  fig.clf()
  ax = seaborn.scatterplot(data = data_frame, x = "Number of RNA sequences", y = "Running time (s)", markers = True)
  fig = ax.get_figure()
  fig.tight_layout()
  image_dir_path = asset_dir_path + "/images"
  fig.savefig(image_dir_path + "/consalifold_running_time_num_rna_seqs.eps", bbox_inches = "tight")
  fig.clf()

def read_consalifold_results(consalifold_output_file_path):
  with open(consalifold_output_file_path) as f:
    lines = f.readlines()
    times = []
    max_seq_lens = []
    nums_rna_seqs = []
    for line in lines:
      split = line.split()
      time = float(split[0])
      max_seq_len = int(split[1])
      num_rna_seqs = int(split[2])
      times.append(time)
      max_seq_lens.append(max_seq_len)
      nums_rna_seqs.append(num_rna_seqs)
    return (times, max_seq_lens, nums_rna_seqs)

def write_consalifold_results(consalifold_results, consalifold_output_file_path):
  with open(consalifold_output_file_path, "w") as f:
    buf = ""
    for consalifold_result in consalifold_results:
      buf += "%f %d %d\n" % (consalifold_result[0], consalifold_result[1], consalifold_result[2])
    f.write(buf)

def bench_consalifold(consalifold_params):
  (sub_thread_num, input_sa_file_path, consalifold_output_dir_path, is_contra_model) = consalifold_params
  consalifold_command = "consalifold %s-g 1 -t " % ("-m contra " if is_contra_model else "") + str(sub_thread_num) + " -i " + input_sa_file_path + " -o " + consalifold_output_dir_path
  begin = time.time()
  utils.run_command(consalifold_command)
  consalifold_elapsed_time = time.time() - begin
  max_seq_len = 0
  sa = AlignIO.read(input_sa_file_path, "clustal")
  for rec in sa:
    seq = str(rec.seq)
    seq_len = len(seq.replace("-", ""))
    if seq_len > max_seq_len:
      max_seq_len = seq_len
  num_rna_seqs = len(sa)
  return (consalifold_elapsed_time, max_seq_len, num_rna_seqs)

if __name__ == "__main__":
  main()
