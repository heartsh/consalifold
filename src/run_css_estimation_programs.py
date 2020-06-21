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
  gammas = [2. ** i for i in range(-7, 11)]
  short_mafft_ginsi_plus_centroidalifold_params_4_bpp_mat = []
  short_mafft_xinsi_plus_centroidalifold_params_4_bpp_mat = []
  short_ref_sa_plus_centroidalifold_params_4_bpp_mat = []
  short_mafft_ginsi_plus_consalifold_params = []
  short_mafft_xinsi_plus_consalifold_params = []
  short_ref_sa_plus_consalifold_params = []
  short_mafft_ginsi_plus_centroidalifold_params = []
  short_mafft_xinsi_plus_centroidalifold_params = []
  short_ref_sa_plus_centroidalifold_params = []
  short_centroidalifold_params_4_elapsed_time = []
  short_mafft_ginsi_plus_rnaalifold_params = []
  short_mafft_xinsi_plus_rnaalifold_params = []
  short_ref_sa_plus_rnaalifold_params = []
  long_mafft_ginsi_plus_centroidalifold_params_4_bpp_mat = []
  long_mafft_xinsi_plus_centroidalifold_params_4_bpp_mat = []
  long_ref_sa_plus_centroidalifold_params_4_bpp_mat = []
  long_mafft_ginsi_plus_consalifold_params = []
  long_mafft_xinsi_plus_consalifold_params = []
  long_ref_sa_plus_consalifold_params = []
  long_mafft_ginsi_plus_centroidalifold_params = []
  long_mafft_xinsi_plus_centroidalifold_params = []
  long_ref_sa_plus_centroidalifold_params = []
  long_centroidalifold_params_4_elapsed_time = []
  long_mafft_ginsi_plus_rnaalifold_params = []
  long_mafft_xinsi_plus_rnaalifold_params = []
  long_ref_sa_plus_rnaalifold_params = []
  short_consalifold_elapsed_time = 0.
  long_consalifold_elapsed_time = 0.
  for data_set in ["short", "long"]:
    # rna_seq_dir_path = asset_dir_path + "/compiled_rna_fams_" + data_set
    rna_seq_dir_path = asset_dir_path + "/compiled_rna_fams_" + data_set + "_4_micro_bench"
    # ref_sa_dir_path = asset_dir_path + "/ref_sas_" + data_set
    ref_sa_dir_path = asset_dir_path + "/ref_sas_" + data_set + "_4_micro_bench"
    mafft_ginsi_dir_path = asset_dir_path + "/mafft_ginsi_" + data_set
    mafft_xinsi_dir_path = asset_dir_path + "/mafft_xinsi_" + data_set
    mafft_ginsi_plus_centroidalifold_bpp_mat_dir_path = asset_dir_path + "/mafft_ginsi_plus_centroidalifold_bpp_mats_" + data_set
    mafft_xinsi_plus_centroidalifold_bpp_mat_dir_path = asset_dir_path + "/mafft_xinsi_plus_centroidalifold_bpp_mats_" + data_set
    ref_sa_plus_centroidalifold_bpp_mat_dir_path = asset_dir_path + "/ref_sa_plus_centroidalifold_bpp_mats_" + data_set
    mafft_ginsi_plus_consalifold_dir_path = asset_dir_path + "/mafft_ginsi_plus_consalifold_" + data_set
    mafft_xinsi_plus_consalifold_dir_path = asset_dir_path + "/mafft_xinsi_plus_consalifold_" + data_set
    ref_sa_plus_consalifold_dir_path = asset_dir_path + "/ref_sa_plus_consalifold_" + data_set
    mafft_ginsi_plus_centroidalifold_dir_path = asset_dir_path + "/mafft_ginsi_plus_centroidalifold_" + data_set
    mafft_xinsi_plus_centroidalifold_dir_path = asset_dir_path + "/mafft_xinsi_plus_centroidalifold_" + data_set
    ref_sa_plus_centroidalifold_dir_path = asset_dir_path + "/ref_sa_plus_centroidalifold_" + data_set
    mafft_ginsi_plus_rnaalifold_dir_path = asset_dir_path + "/mafft_ginsi_plus_rnaalifold_" + data_set
    mafft_xinsi_plus_rnaalifold_dir_path = asset_dir_path + "/mafft_xinsi_plus_rnaalifold_" + data_set
    ref_sa_plus_rnaalifold_dir_path = asset_dir_path + "/ref_sa_plus_rnaalifold_" + data_set
    if not os.path.isdir(mafft_ginsi_dir_path):
      os.mkdir(mafft_ginsi_dir_path)
    if not os.path.isdir(mafft_xinsi_dir_path):
      os.mkdir(mafft_xinsi_dir_path)
    if not os.path.isdir(ref_sa_dir_path):
      os.mkdir(ref_sa_dir_path)
    if not os.path.isdir(mafft_ginsi_plus_centroidalifold_bpp_mat_dir_path):
      os.mkdir(mafft_ginsi_plus_centroidalifold_bpp_mat_dir_path)
    if not os.path.isdir(mafft_xinsi_plus_centroidalifold_bpp_mat_dir_path):
      os.mkdir(mafft_xinsi_plus_centroidalifold_bpp_mat_dir_path)
    if not os.path.isdir(ref_sa_plus_centroidalifold_bpp_mat_dir_path):
      os.mkdir(ref_sa_plus_centroidalifold_bpp_mat_dir_path)
    if not os.path.isdir(mafft_ginsi_plus_consalifold_dir_path):
      os.mkdir(mafft_ginsi_plus_consalifold_dir_path)
    if not os.path.isdir(mafft_xinsi_plus_consalifold_dir_path):
      os.mkdir(mafft_xinsi_plus_consalifold_dir_path)
    if not os.path.isdir(ref_sa_plus_consalifold_dir_path):
      os.mkdir(ref_sa_plus_consalifold_dir_path)
    if not os.path.isdir(mafft_ginsi_plus_centroidalifold_dir_path):
      os.mkdir(mafft_ginsi_plus_centroidalifold_dir_path)
    if not os.path.isdir(mafft_xinsi_plus_centroidalifold_dir_path):
      os.mkdir(mafft_xinsi_plus_centroidalifold_dir_path)
    if not os.path.isdir(ref_sa_plus_centroidalifold_dir_path):
      os.mkdir(ref_sa_plus_centroidalifold_dir_path)
    if not os.path.isdir(mafft_ginsi_plus_rnaalifold_dir_path):
      os.mkdir(mafft_ginsi_plus_rnaalifold_dir_path)
    if not os.path.isdir(mafft_xinsi_plus_rnaalifold_dir_path):
      os.mkdir(mafft_xinsi_plus_rnaalifold_dir_path)
    if not os.path.isdir(ref_sa_plus_rnaalifold_dir_path):
      os.mkdir(ref_sa_plus_rnaalifold_dir_path)
    for rna_seq_file in os.listdir(rna_seq_dir_path):
      if not rna_seq_file.endswith(".fa"):
        continue
      rna_seq_file_path = os.path.join(rna_seq_dir_path, rna_seq_file)
      (rna_family_name, extension) = os.path.splitext(rna_seq_file)
      mafft_ginsi_plus_centroidalifold_bpp_mat_file_path = os.path.join(mafft_ginsi_plus_centroidalifold_bpp_mat_dir_path, rna_family_name + ".dat")
      mafft_xinsi_plus_centroidalifold_bpp_mat_file_path = os.path.join(mafft_xinsi_plus_centroidalifold_bpp_mat_dir_path, rna_family_name + ".dat")
      ref_sa_plus_centroidalifold_bpp_mat_file_path = os.path.join(ref_sa_plus_centroidalifold_bpp_mat_dir_path, rna_family_name + ".dat")
      mafft_ginsi_output_file_path = os.path.join(mafft_ginsi_dir_path, rna_family_name + ".aln")
      mafft_xinsi_output_file_path = os.path.join(mafft_xinsi_dir_path, rna_family_name + ".aln")
      run_mafft((rna_seq_file_path, mafft_ginsi_output_file_path, "ginsi", num_of_threads))
      run_mafft((rna_seq_file_path, mafft_xinsi_output_file_path, "xinsi", num_of_threads))
      ref_sa_file_path = os.path.join(ref_sa_dir_path, rna_family_name + ".aln")
      mafft_ginsi_plus_centroidalifold_command = "centroid_alifold -e McCaskill -w 0 -e Alifold -w 1 --posteriors 0 --posteriors-output " + mafft_ginsi_plus_centroidalifold_bpp_mat_file_path + " " + mafft_ginsi_output_file_path
      mafft_xinsi_plus_centroidalifold_command = "centroid_alifold -e McCaskill -w 0 -e Alifold -w 1 --posteriors 0 --posteriors-output " + mafft_xinsi_plus_centroidalifold_bpp_mat_file_path + " " + mafft_xinsi_output_file_path
      ref_sa_plus_centroidalifold_command = "centroid_alifold -e McCaskill -w 0 -e Alifold -w 1 --posteriors 0 --posteriors-output " + ref_sa_plus_centroidalifold_bpp_mat_file_path + " " + ref_sa_file_path
      if data_set == "short":
        short_mafft_ginsi_plus_centroidalifold_params_4_bpp_mat.insert(0, mafft_ginsi_plus_centroidalifold_command)
        short_mafft_xinsi_plus_centroidalifold_params_4_bpp_mat.insert(0, mafft_xinsi_plus_centroidalifold_command)
        short_ref_sa_plus_centroidalifold_params_4_bpp_mat.insert(0, ref_sa_plus_centroidalifold_command)
      else:
        long_mafft_ginsi_plus_centroidalifold_params_4_bpp_mat.insert(0, mafft_ginsi_plus_centroidalifold_command)
        long_mafft_xinsi_plus_centroidalifold_params_4_bpp_mat.insert(0, mafft_xinsi_plus_centroidalifold_command)
        long_ref_sa_plus_centroidalifold_params_4_bpp_mat.insert(0, ref_sa_plus_centroidalifold_command)
  pool = multiprocessing.Pool(num_of_threads)
  pool.map(utils.run_command, short_mafft_ginsi_plus_centroidalifold_params_4_bpp_mat)
  begin = time.time()
  pool.map(utils.run_command, short_mafft_xinsi_plus_centroidalifold_params_4_bpp_mat)
  short_consalifold_elapsed_time += time.time() - begin
  pool.map(utils.run_command, short_ref_sa_plus_centroidalifold_params_4_bpp_mat)
  pool.map(utils.run_command, long_mafft_ginsi_plus_centroidalifold_params_4_bpp_mat)
  begin = time.time()
  pool.map(utils.run_command, long_mafft_xinsi_plus_centroidalifold_params_4_bpp_mat)
  long_consalifold_elapsed_time += time.time() - begin
  pool.map(utils.run_command, long_ref_sa_plus_centroidalifold_params_4_bpp_mat)
  sub_thread_num = 4
  for data_set in ["short", "long"]:
    # rna_seq_dir_path = asset_dir_path + "/compiled_rna_fams_" + data_set
    rna_seq_dir_path = asset_dir_path + "/compiled_rna_fams_" + data_set + "_4_micro_bench"
    # ref_sa_dir_path = asset_dir_path + "/ref_sas_" + data_set
    ref_sa_dir_path = asset_dir_path + "/ref_sas_" + data_set + "_4_micro_bench"
    mafft_ginsi_dir_path = asset_dir_path + "/mafft_ginsi_" + data_set
    mafft_xinsi_dir_path = asset_dir_path + "/mafft_xinsi_" + data_set
    mafft_ginsi_plus_centroidalifold_bpp_mat_dir_path = asset_dir_path + "/mafft_ginsi_plus_centroidalifold_bpp_mats_" + data_set
    mafft_xinsi_plus_centroidalifold_bpp_mat_dir_path = asset_dir_path + "/mafft_xinsi_plus_centroidalifold_bpp_mats_" + data_set
    ref_sa_plus_centroidalifold_bpp_mat_dir_path = asset_dir_path + "/ref_sa_plus_centroidalifold_bpp_mats_" + data_set
    mafft_ginsi_plus_consalifold_dir_path = asset_dir_path + "/mafft_ginsi_plus_consalifold_" + data_set
    mafft_xinsi_plus_consalifold_dir_path = asset_dir_path + "/mafft_xinsi_plus_consalifold_" + data_set
    ref_sa_plus_consalifold_dir_path = asset_dir_path + "/ref_sa_plus_consalifold_" + data_set
    mafft_ginsi_plus_centroidalifold_dir_path = asset_dir_path + "/mafft_ginsi_plus_centroidalifold_" + data_set
    mafft_xinsi_plus_centroidalifold_dir_path = asset_dir_path + "/mafft_xinsi_plus_centroidalifold_" + data_set
    ref_sa_plus_centroidalifold_dir_path = asset_dir_path + "/ref_sa_plus_centroidalifold_" + data_set
    mafft_ginsi_plus_rnaalifold_dir_path = asset_dir_path + "/mafft_ginsi_plus_rnaalifold_" + data_set
    mafft_xinsi_plus_rnaalifold_dir_path = asset_dir_path + "/mafft_xinsi_plus_rnaalifold_" + data_set
    ref_sa_plus_rnaalifold_dir_path = asset_dir_path + "/ref_sa_plus_rnaalifold_" + data_set
    for rna_seq_file in os.listdir(rna_seq_dir_path):
      if not rna_seq_file.endswith(".fa"):
        continue
      rna_seq_file_path = os.path.join(rna_seq_dir_path, rna_seq_file)
      (rna_family_name, extension) = os.path.splitext(rna_seq_file)
      mafft_ginsi_plus_centroidalifold_bpp_mat_file_path = os.path.join(mafft_ginsi_plus_centroidalifold_bpp_mat_dir_path, rna_family_name + ".dat")
      mafft_xinsi_plus_centroidalifold_bpp_mat_file_path = os.path.join(mafft_xinsi_plus_centroidalifold_bpp_mat_dir_path, rna_family_name + ".dat")
      ref_sa_plus_centroidalifold_bpp_mat_file_path = os.path.join(ref_sa_plus_centroidalifold_bpp_mat_dir_path, rna_family_name + ".dat")
      mafft_ginsi_output_file_path = os.path.join(mafft_ginsi_dir_path, rna_family_name + ".aln")
      mafft_xinsi_output_file_path = os.path.join(mafft_xinsi_dir_path, rna_family_name + ".aln")
      ref_sa_file_path = os.path.join(ref_sa_dir_path, rna_family_name + ".aln")
      mafft_ginsi_plus_consalifold_output_dir_path = os.path.join(mafft_ginsi_plus_consalifold_dir_path, rna_family_name)
      mafft_xinsi_plus_consalifold_output_dir_path = os.path.join(mafft_xinsi_plus_consalifold_dir_path, rna_family_name)
      ref_sa_plus_consalifold_output_dir_path = os.path.join(ref_sa_plus_consalifold_dir_path, rna_family_name)
      mafft_ginsi_plus_centroidalifold_output_dir_path = os.path.join(mafft_ginsi_plus_centroidalifold_dir_path, rna_family_name)
      mafft_xinsi_plus_centroidalifold_output_dir_path = os.path.join(mafft_xinsi_plus_centroidalifold_dir_path, rna_family_name)
      ref_sa_plus_centroidalifold_output_dir_path = os.path.join(ref_sa_plus_centroidalifold_dir_path, rna_family_name)
      if not os.path.isdir(mafft_ginsi_plus_consalifold_output_dir_path):
        os.mkdir(mafft_ginsi_plus_consalifold_output_dir_path)
      if not os.path.isdir(mafft_xinsi_plus_consalifold_output_dir_path):
        os.mkdir(mafft_xinsi_plus_consalifold_output_dir_path)
      if not os.path.isdir(ref_sa_plus_consalifold_output_dir_path):
        os.mkdir(ref_sa_plus_consalifold_output_dir_path)
      if not os.path.isdir(mafft_ginsi_plus_centroidalifold_output_dir_path):
        os.mkdir(mafft_ginsi_plus_centroidalifold_output_dir_path)
      if not os.path.isdir(mafft_xinsi_plus_centroidalifold_output_dir_path):
        os.mkdir(mafft_xinsi_plus_centroidalifold_output_dir_path)
      if not os.path.isdir(ref_sa_plus_centroidalifold_output_dir_path):
        os.mkdir(ref_sa_plus_centroidalifold_output_dir_path)
      mafft_ginsi_plus_consalifold_command = "consalifold -t " + str(sub_thread_num) + " -i " + rna_seq_file_path + " -a " + mafft_ginsi_output_file_path + " -c " + mafft_ginsi_plus_centroidalifold_bpp_mat_file_path + " -o " + mafft_ginsi_plus_consalifold_output_dir_path
      mafft_xinsi_plus_consalifold_command = "consalifold -t " + str(sub_thread_num) + " -i " + rna_seq_file_path + " -a " + mafft_xinsi_output_file_path + " -c " + mafft_xinsi_plus_centroidalifold_bpp_mat_file_path + " -o " + mafft_xinsi_plus_consalifold_output_dir_path
      ref_sa_plus_consalifold_command = "consalifold -t " + str(sub_thread_num) + " -i " + rna_seq_file_path + " -a " + ref_sa_file_path + " -c " + ref_sa_plus_centroidalifold_bpp_mat_file_path + " -o " + ref_sa_plus_consalifold_output_dir_path
      if data_set == "short":
        short_mafft_ginsi_plus_consalifold_params.insert(0, mafft_ginsi_plus_consalifold_command)
        short_mafft_xinsi_plus_consalifold_params.insert(0, mafft_xinsi_plus_consalifold_command)
        short_ref_sa_plus_consalifold_params.insert(0, ref_sa_plus_consalifold_command)
      else:
        long_mafft_ginsi_plus_consalifold_params.insert(0, mafft_ginsi_plus_consalifold_command)
        long_mafft_xinsi_plus_consalifold_params.insert(0, mafft_xinsi_plus_consalifold_command)
        long_ref_sa_plus_consalifold_params.insert(0, ref_sa_plus_consalifold_command)
      output_file = rna_family_name + ".sth"
      mafft_ginsi_plus_rnaalifold_output_file_path = os.path.join(mafft_ginsi_plus_rnaalifold_dir_path, output_file)
      mafft_xinsi_plus_rnaalifold_output_file_path = os.path.join(mafft_xinsi_plus_rnaalifold_dir_path, output_file)
      ref_sa_plus_rnaalifold_output_file_path = os.path.join(ref_sa_plus_rnaalifold_dir_path, output_file)
      if data_set == "short":
        short_mafft_ginsi_plus_rnaalifold_params.insert(0, (mafft_ginsi_output_file_path, mafft_ginsi_plus_rnaalifold_output_file_path))
        short_mafft_xinsi_plus_rnaalifold_params.insert(0, (mafft_xinsi_output_file_path, mafft_xinsi_plus_rnaalifold_output_file_path))
        short_ref_sa_plus_rnaalifold_params.insert(0, (ref_sa_file_path, ref_sa_plus_rnaalifold_output_file_path))
      else:
        long_mafft_ginsi_plus_rnaalifold_params.insert(0, (mafft_ginsi_output_file_path, mafft_ginsi_plus_rnaalifold_output_file_path))
        long_mafft_xinsi_plus_rnaalifold_params.insert(0, (mafft_xinsi_output_file_path, mafft_xinsi_plus_rnaalifold_output_file_path))
        long_ref_sa_plus_rnaalifold_params.insert(0, (ref_sa_file_path, ref_sa_plus_rnaalifold_output_file_path))
      for gamma in gammas:
        gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
        output_file = "gamma=" + gamma_str + ".sth"
        mafft_ginsi_plus_centroidalifold_output_file_path = os.path.join(mafft_ginsi_plus_centroidalifold_output_dir_path, output_file)
        mafft_xinsi_plus_centroidalifold_output_file_path = os.path.join(mafft_xinsi_plus_centroidalifold_output_dir_path, output_file)
        ref_sa_plus_centroidalifold_output_file_path = os.path.join(ref_sa_plus_centroidalifold_output_dir_path, output_file)
        if data_set == "short":
          short_mafft_ginsi_plus_centroidalifold_params.insert(0, (mafft_ginsi_output_file_path, mafft_ginsi_plus_centroidalifold_output_file_path, gamma_str))
          short_mafft_xinsi_plus_centroidalifold_params.insert(0, (mafft_xinsi_output_file_path, mafft_xinsi_plus_centroidalifold_output_file_path, gamma_str))
          short_ref_sa_plus_centroidalifold_params.insert(0, (ref_sa_file_path, ref_sa_plus_centroidalifold_output_file_path, gamma_str))
          if gamma == 1:
            short_centroidalifold_params_4_elapsed_time.insert(0, (mafft_xinsi_output_file_path, mafft_xinsi_plus_centroidalifold_output_file_path, gamma_str))
        else:
          long_mafft_ginsi_plus_centroidalifold_params.insert(0, (mafft_ginsi_output_file_path, mafft_ginsi_plus_centroidalifold_output_file_path, gamma_str))
          long_mafft_xinsi_plus_centroidalifold_params.insert(0, (mafft_xinsi_output_file_path, mafft_xinsi_plus_centroidalifold_output_file_path, gamma_str))
          long_ref_sa_plus_centroidalifold_params.insert(0, (ref_sa_file_path, ref_sa_plus_centroidalifold_output_file_path, gamma_str))
          if gamma == 1:
            long_centroidalifold_params_4_elapsed_time.insert(0, (mafft_xinsi_output_file_path, mafft_xinsi_plus_centroidalifold_output_file_path, gamma_str))
  pool = multiprocessing.Pool(int(num_of_threads / sub_thread_num))
  pool.map(utils.run_command, short_mafft_ginsi_plus_consalifold_params)
  begin = time.time()
  pool.map(utils.run_command, short_mafft_xinsi_plus_consalifold_params)
  short_consalifold_elapsed_time += time.time() - begin
  pool.map(utils.run_command, short_ref_sa_plus_consalifold_params)
  pool.map(utils.run_command, long_mafft_ginsi_plus_consalifold_params)
  begin = time.time()
  pool.map(utils.run_command, long_mafft_xinsi_plus_consalifold_params)
  long_consalifold_elapsed_time += time.time() - begin
  pool.map(utils.run_command, long_ref_sa_plus_consalifold_params)
  pool = multiprocessing.Pool(num_of_threads)
  pool.map(run_centroidalifold, short_mafft_ginsi_plus_centroidalifold_params)
  pool.map(run_centroidalifold, short_mafft_xinsi_plus_centroidalifold_params)
  pool.map(run_centroidalifold, short_ref_sa_plus_centroidalifold_params)
  begin = time.time()
  pool.map(run_centroidalifold, short_centroidalifold_params_4_elapsed_time)
  short_centroidalifold_elapsed_time = time.time() - begin
  pool.map(run_rnaalifold, short_mafft_ginsi_plus_rnaalifold_params)
  begin = time.time()
  pool.map(run_rnaalifold, short_mafft_xinsi_plus_rnaalifold_params)
  short_rnaalifold_elapsed_time = time.time() - begin
  pool.map(run_rnaalifold, short_ref_sa_plus_rnaalifold_params)
  pool.map(run_centroidalifold, long_mafft_ginsi_plus_centroidalifold_params)
  pool.map(run_centroidalifold, long_mafft_xinsi_plus_centroidalifold_params)
  pool.map(run_centroidalifold, long_ref_sa_plus_centroidalifold_params)
  begin = time.time()
  pool.map(run_centroidalifold, long_centroidalifold_params_4_elapsed_time)
  long_centroidalifold_elapsed_time = time.time() - begin
  pool.map(run_rnaalifold, long_mafft_ginsi_plus_rnaalifold_params)
  begin = time.time()
  pool.map(run_rnaalifold, long_mafft_xinsi_plus_rnaalifold_params)
  long_rnaalifold_elapsed_time = time.time() - begin
  pool.map(run_rnaalifold, long_ref_sa_plus_rnaalifold_params)
  print("The elapsed time of the ConsAlifold program for test set \"short\" = %f [s]." % short_consalifold_elapsed_time)
  print("The elapsed time of the CentroidAlifold program for test set \"short\" = %f [s]." % short_centroidalifold_elapsed_time)
  print("The elapsed time of the RNAalifold program for test set \"short\" = %f [s]." % short_rnaalifold_elapsed_time)
  print("The elapsed time of the ConsAlifold program for test set \"long\" = %f [s]." % long_consalifold_elapsed_time)
  print("The elapsed time of the CentroidAlifold program for test set \"long\" = %f [s]." % long_centroidalifold_elapsed_time)
  print("The elapsed time of the RNAalifold program for test set \"long\" = %f [s]." % long_rnaalifold_elapsed_time)

def run_mafft(mafft_params):
  (rna_seq_file_path, mafft_output_file_path, mafft_type, num_of_threads) = mafft_params
  mafft_command = "mafft-%s --thread %d --quiet " % (mafft_type, num_of_threads) + "--clustalout " + rna_seq_file_path + " > " + mafft_output_file_path
  utils.run_command(mafft_command)

def run_centroidalifold(centroidalifold_params):
  (sa_file_path, centroidalifold_output_file_path, gamma_str) = centroidalifold_params
  centroidalifold_command = "centroid_alifold " + sa_file_path + " -g " + gamma_str
  (output, _, _) = utils.run_command(centroidalifold_command)
  css = str(output).strip().split("\\n")[2].split()[0]
  sta = AlignIO.read(sa_file_path, "clustal")
  AlignIO.write(sta, centroidalifold_output_file_path, "stockholm")
  sta = AlignIO.read(centroidalifold_output_file_path, "stockholm")
  sta.column_annotations["secondary_structure"] = css
  AlignIO.write(sta, centroidalifold_output_file_path, "stockholm")

def run_rnaalifold(rnaalifold_params):
  (sa_file_path, rnaalifold_output_file_path) = rnaalifold_params
  rnaalifold_command = "RNAalifold -q --noPS " + sa_file_path
  (output, _, _) = utils.run_command(rnaalifold_command)
  css = str(output).strip().split("\\n")[1].split()[0]
  sta = AlignIO.read(sa_file_path, "clustal")
  AlignIO.write(sta, rnaalifold_output_file_path, "stockholm")
  sta = AlignIO.read(rnaalifold_output_file_path, "stockholm")
  sta.column_annotations["secondary_structure"] = css
  AlignIO.write(sta, rnaalifold_output_file_path, "stockholm")

if __name__ == "__main__":
  main()
