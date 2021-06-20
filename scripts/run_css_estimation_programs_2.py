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

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  num_of_threads = multiprocessing.cpu_count()
  temp_dir_path = "/tmp/run_css_estimation_programs_2_%s" % datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d_%H:%M:%S')
  if not os.path.isdir(temp_dir_path):
    os.mkdir(temp_dir_path)
  gammas = [2. ** i for i in range(-7, 11)]
  centroidhomfold_params = []
  rnaspa_params = []
  locarna_params = []
  raf_params = []
  turbofold_params = []
  centroidhomfold_params_4_running_time = []
  turbofold_params_4_running_time = []
  rna_seq_dir_path = asset_dir_path + "/compiled_rna_fams_test"
  centroidhomfold_dir_path = asset_dir_path + "/centroidhomfold"
  rnaspa_dir_path = asset_dir_path + "/rnaspa"
  locarna_dir_path = asset_dir_path + "/locarna"
  raf_dir_path = asset_dir_path + "/raf"
  turbofold_dir_path = asset_dir_path + "/turbofold"
  if not os.path.isdir(centroidhomfold_dir_path):
    os.mkdir(centroidhomfold_dir_path)
  if not os.path.isdir(rnaspa_dir_path):
    os.mkdir(rnaspa_dir_path)
  if not os.path.isdir(locarna_dir_path):
    os.mkdir(locarna_dir_path)
  if not os.path.isdir(raf_dir_path):
    os.mkdir(raf_dir_path)
  if not os.path.isdir(turbofold_dir_path):
    os.mkdir(turbofold_dir_path)
  sub_thread_num = 4
  for rna_seq_file in os.listdir(rna_seq_dir_path):
    if not rna_seq_file.endswith(".fa"):
      continue
    rna_seq_file_path = os.path.join(rna_seq_dir_path, rna_seq_file)
    (rna_family_name, extension) = os.path.splitext(rna_seq_file)
    centroidhomfold_output_dir_path = os.path.join(centroidhomfold_dir_path, rna_family_name)
    rnaspa_output_file_path = os.path.join(rnaspa_dir_path, rna_family_name + ".fa")
    locarna_output_dir_path = os.path.join(locarna_dir_path, rna_family_name)
    raf_output_file_path = os.path.join(raf_dir_path, rna_family_name + ".sth")
    turbofold_output_dir_path = os.path.join(turbofold_dir_path, rna_family_name)
    if not os.path.isdir(centroidhomfold_output_dir_path):
      os.mkdir(centroidhomfold_output_dir_path)
    if not os.path.isdir(locarna_output_dir_path):
      os.mkdir(locarna_output_dir_path)
    if not os.path.isdir(turbofold_output_dir_path):
      os.mkdir(turbofold_output_dir_path)
    rnaspa_params.insert(0, (rna_seq_file_path, rnaspa_output_file_path))
    locarna_params.insert(0, (rna_seq_file_path, locarna_output_dir_path))
    raf_params.insert(0, (rna_seq_file_path, raf_output_file_path))
    for gamma in gammas:
      gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
      output_file = "gamma=" + gamma_str + ".sth"
      output_file_2 = "gamma=" + gamma_str + ".fa"
      centroidhomfold_output_file_path = os.path.join(centroidhomfold_output_dir_path, output_file_2)
      centroidhomfold_params.insert(0, (rna_seq_file_path, centroidhomfold_output_file_path, gamma_str, temp_dir_path))
      if gamma == 1:
        centroidhomfold_params_4_running_time.insert(0, (rna_seq_file_path, centroidhomfold_output_file_path, gamma_str, temp_dir_path))
      turbofold_output_file_path = os.path.join(turbofold_output_dir_path, output_file_2)
      turbofold_params.insert(0, (rna_seq_file_path, turbofold_output_file_path, gamma, temp_dir_path, rna_family_name))
      if gamma == 1:
        turbofold_params_4_running_time.insert(0, (rna_seq_file_path, turbofold_output_file_path, gamma, temp_dir_path, rna_family_name))
  pool = multiprocessing.Pool(num_of_threads)
  pool.map(run_centroidhomfold, centroidhomfold_params)
  begin = time.time()
  pool.map(run_rnaspa, rnaspa_params)
  rnaspa_elapsed_time = time.time() - begin
  begin = time.time()
  pool.map(run_locarna, locarna_params)
  locarna_elapsed_time = time.time() - begin
  begin = time.time()
  pool.map(run_raf, raf_params)
  raf_elapsed_time = time.time() - begin
  pool.map(run_turbofold, turbofold_params)
  begin = time.time()
  pool.map(run_centroidhomfold, centroidhomfold_params_4_running_time)
  centroidhomfold_elapsed_time = time.time() - begin
  begin = time.time()
  pool.map(run_turbofold, turbofold_params_4_running_time)
  turbofold_elapsed_time = time.time() - begin
  print("The elapsed time of CentroidHomfold = %f [s]." % centroidhomfold_elapsed_time)
  print("The elapsed time of RNAspa = %f [s]." % rnaspa_elapsed_time)
  print("The elapsed time of LocARNA = %f [s]." % locarna_elapsed_time)
  print("The elapsed time of RAF = %f [s]." % raf_elapsed_time)
  print("The elapsed time of TurboFold = %f [s]." % turbofold_elapsed_time)

def run_locarna(locarna_params):
  (rna_file_path, locarna_output_dir_path) = locarna_params
  locarna_command = "mlocarna %s --keep-sequence-order --tgtdir %s --stockholm --consensus-structure alifold" % (rna_file_path, locarna_output_dir_path)
  utils.run_command(locarna_command)
  sa_file_path = locarna_output_dir_path + "/results/result.aln"
  sa = AlignIO.read(sa_file_path, "clustal")
  sta_file_path = locarna_output_dir_path + "/results/result.stk"
  sta = AlignIO.read(sta_file_path, "stockholm")
  css_string = sta.column_annotations["secondary_structure"]
  sa.column_annotations["secondary_structure"] = css_string
  AlignIO.write(sa, sta_file_path, "stockholm")

def run_turbofold(turbofold_params):
  (rna_file_path, turbofold_output_file_path, gamma, temp_dir_path, rna_family_name) = turbofold_params
  recs = [rec for rec in SeqIO.parse(rna_file_path, "fasta")]
  rec_lens = [len(rec) for rec in recs]
  rec_seq_len = len(recs)
  turbofold_temp_dir_path = "%s/%s_gamma=%d" % (temp_dir_path, rna_family_name, gamma)
  if not os.path.isdir(turbofold_temp_dir_path):
    os.mkdir(turbofold_temp_dir_path)
  turbofold_config_file_contents = "InSeq = {"
  for i in range(rec_seq_len):
    turbofold_config_file_contents += "%s/%d.fasta;" % (turbofold_temp_dir_path, i)
  turbofold_config_file_contents += "}\nOutCT = {"
  for i in range(rec_seq_len):
    turbofold_config_file_contents += "%s/%d.ct;" % (turbofold_temp_dir_path, i)
  turbofold_config_file_contents += "}\nIterations = 3\nMode = MEA\nMeaGamma = %f" % gamma
  turbofold_config_file_path = os.path.join(turbofold_temp_dir_path, "turbofold_config.dat")
  turbofold_config_file = open(turbofold_config_file_path, "w")
  turbofold_config_file.write(turbofold_config_file_contents)
  turbofold_config_file.close()
  for (i, rec) in enumerate(recs):
    SeqIO.write([rec], open(os.path.join(turbofold_temp_dir_path, "%d.fasta" % i), "w"), "fasta")
  turbofold_command = "TurboFold " + turbofold_config_file_path
  utils.run_command(turbofold_command)
  turbofold_output_file_contents = ""
  all_files_exist = True
  for i in range(rec_seq_len):
    ct_file_path = os.path.join(turbofold_temp_dir_path, "%d.ct" % i)
    if os.path.exists(ct_file_path):
      ss_string = read_ct_file(ct_file_path)
      turbofold_output_file_contents += ">%d\n%s\n\n" % (i, ss_string)
    else:
      all_files_exist = False
      turbofold_output_file_contents = ""
      break
  if not all_files_exist:
    print("Some output files are empty. TurboFold is retried with # iterations = 1.")
    turbofold_config_file_contents = "InSeq = {"
    for i in range(rec_seq_len):
      turbofold_config_file_contents += "%s/%d.fasta;" % (turbofold_temp_dir_path, i)
    turbofold_config_file_contents += "}\nOutCT = {"
    for i in range(rec_seq_len):
      turbofold_config_file_contents += "%s/%d.ct;" % (turbofold_temp_dir_path, i)
    turbofold_config_file_contents += "}\nIterations = 1\nMode = MEA\nMeaGamma = %f" % gamma
    turbofold_config_file_path = os.path.join(turbofold_temp_dir_path, "turbofold_config.dat")
    turbofold_config_file = open(turbofold_config_file_path, "w")
    turbofold_config_file.write(turbofold_config_file_contents)
    turbofold_config_file.close()
    utils.run_command(turbofold_command)
    for i in range(rec_seq_len):
      ct_file_path = os.path.join(turbofold_temp_dir_path, "%d.ct" % i)
      ss_string = read_ct_file(ct_file_path)
      turbofold_output_file_contents += ">%d\n%s\n\n" % (i, ss_string)
  turbofold_output_file = open(turbofold_output_file_path, "w")
  turbofold_output_file.write(turbofold_output_file_contents)
  turbofold_output_file.close()

def read_ct_file(ct_file_path):
  ct_file = open(ct_file_path, "r")
  lines = ct_file.readlines()
  seq_len = int(lines[0].split()[0])
  ss_string = ["." for i in range(seq_len)]
  num_of_lines = len(lines)
  for line in lines[1 : num_of_lines]:
    if "ENERGY" in line:
      break
    substrings = line.split()
    index_1 = int(substrings[0])
    index_2 = int(substrings[4])
    if index_2 == 0 or index_1 >= index_2:
      continue
    ss_string[index_1 - 1] = "("
    ss_string[index_2 - 1] = ")"
  return "".join(ss_string)

def run_centroidhomfold(centroidhomfold_params):
  (rna_file_path, centroidhomfold_output_file_path, gamma_str, temp_dir_path) = centroidhomfold_params
  recs = [rec for rec in SeqIO.parse(rna_file_path, "fasta")]
  rec_seq_len = len(recs)
  centroidhomfold_output_file = open(centroidhomfold_output_file_path, "w+")
  centroidhomfold_output_buf = ""
  basename = os.path.basename(rna_file_path)
  (rna_family_name, extension) = os.path.splitext(basename)
  seq_file_path = os.path.join(temp_dir_path, "seqs_4_%s_and_gamma=%s.fa" % (rna_family_name, gamma_str))
  hom_seq_file_path = os.path.join(temp_dir_path, "hom_seqs_4_%s_and_gamma=%s.fa" % (rna_family_name, gamma_str))
  for (i, rec) in enumerate(recs):
    SeqIO.write([rec], open(seq_file_path, "w"), "fasta")
    hom_recs = [rec for (j, rec) in enumerate(recs) if j != i]
    SeqIO.write(recs, open(hom_seq_file_path, "w"), "fasta")
    centroidhomfold_command = "centroid_homfold " + seq_file_path + " -H " + hom_seq_file_path + " -g " + gamma_str
    (output, _, _) = utils.run_command(centroidhomfold_command)
    centroidhomfold_output_buf += ">%d\n%s\n\n" % (i, str(output).split("\\n")[2].split()[0])
  centroidhomfold_output_file.write(centroidhomfold_output_buf)
  centroidhomfold_output_file.close()

def run_rnaspa(rnaspa_params):
  (rna_file_path, rnaspa_output_file_path) = rnaspa_params
  rnaspa_command = "RNAspa %s" % rna_file_path
  (output, _, _) = utils.run_command(rnaspa_command)
  sss = [line.strip() for line in str(output).split("\\n") if line.startswith(".") or line.startswith("(")]
  rnaspa_output_file = open(rnaspa_output_file_path, "w+")
  rnaspa_output_buf = ""
  for (i, ss) in enumerate(sss):
    rnaspa_output_buf += ">%d\n%s\n\n" % (i, ss)
  rnaspa_output_file.write(rnaspa_output_buf)
  rnaspa_output_file.close()

def run_raf(raf_params):
  (rna_file_path, raf_output_file_path) = raf_params
  raf_command = "raf predict " + rna_file_path
  (output, _, _) = utils.run_command(raf_command)
  raf_output_file = open(raf_output_file_path, "w+")
  raf_output_file.write(output.decode())
  raf_output_file.close()
  sta = AlignIO.read(raf_output_file_path, "fasta")
  recs = sta[:-1]
  new_sta = AlignIO.MultipleSeqAlignment(recs)
  new_sta.column_annotations["secondary_structure"] = str(sta[-1].seq)
  AlignIO.write(new_sta, raf_output_file_path, "stockholm")

if __name__ == "__main__":
  main()
