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
  mafft_ginsi_dir_path = asset_dir_path + "/mafft_ginsi"
  mafft_xinsi_dir_path = asset_dir_path + "/mafft_xinsi"
  ref_sa_dir_path = asset_dir_path + "/ref_sas"
  mafft_ginsi_plus_centroidalifold_bpp_mat_dir_path = asset_dir_path + "/mafft_ginsi_plus_centroidalifold_bpp_mats"
  mafft_xinsi_plus_centroidalifold_bpp_mat_dir_path = asset_dir_path + "/mafft_xinsi_plus_centroidalifold_bpp_mats"
  ref_sa_plus_centroidalifold_bpp_mat_dir_path = asset_dir_path + "/ref_sa_plus_centroidalifold_bpp_mats"
  mafft_ginsi_plus_phyloalifold_dir_path = asset_dir_path + "/mafft_ginsi_plus_phyloalifold"
  mafft_xinsi_plus_phyloalifold_dir_path = asset_dir_path + "/mafft_xinsi_plus_phyloalifold"
  ref_sa_plus_phyloalifold_dir_path = asset_dir_path + "/ref_sa_plus_phyloalifold"
  mafft_ginsi_plus_centroidalifold_dir_path = asset_dir_path + "/mafft_ginsi_plus_centroidalifold"
  mafft_xinsi_plus_centroidalifold_dir_path = asset_dir_path + "/mafft_xinsi_plus_centroidalifold"
  ref_sa_plus_centroidalifold_dir_path = asset_dir_path + "/ref_sa_plus_centroidalifold"
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
  if not os.path.isdir(mafft_ginsi_plus_phyloalifold_dir_path):
    os.mkdir(mafft_ginsi_plus_phyloalifold_dir_path)
  if not os.path.isdir(mafft_xinsi_plus_phyloalifold_dir_path):
    os.mkdir(mafft_xinsi_plus_phyloalifold_dir_path)
  if not os.path.isdir(ref_sa_plus_phyloalifold_dir_path):
    os.mkdir(ref_sa_plus_phyloalifold_dir_path)
  if not os.path.isdir(mafft_ginsi_plus_centroidalifold_dir_path):
    os.mkdir(mafft_ginsi_plus_centroidalifold_dir_path)
  if not os.path.isdir(mafft_xinsi_plus_centroidalifold_dir_path):
    os.mkdir(mafft_xinsi_plus_centroidalifold_dir_path)
  if not os.path.isdir(ref_sa_plus_centroidalifold_dir_path):
    os.mkdir(ref_sa_plus_centroidalifold_dir_path)
  rna_seq_dir_path = asset_dir_path + "/sampled_rna_fams"
  gammas = [2. ** i for i in range(-7, 11)]
  mafft_ginsi_plus_centroidalifold_params_4_bpp_mat = []
  mafft_xinsi_plus_centroidalifold_params_4_bpp_mat = []
  ref_sa_plus_centroidalifold_params_4_bpp_mat = []
  mafft_ginsi_plus_phyloalifold_params = []
  mafft_xinsi_plus_phyloalifold_params = []
  ref_sa_plus_phyloalifold_params = []
  mafft_ginsi_plus_centroidalifold_params = []
  mafft_xinsi_plus_centroidalifold_params = []
  ref_sa_plus_centroidalifold_params = []
  centroidalifold_params_4_elapsed_time = []
  phyloalifold_elapsed_time = 0.
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
    mafft_ginsi_plus_centroidalifold_params_4_bpp_mat.insert(0, mafft_ginsi_plus_centroidalifold_command)
    mafft_xinsi_plus_centroidalifold_command = "centroid_alifold -e McCaskill -w 0 -e Alifold -w 1 --posteriors 0 --posteriors-output " + mafft_xinsi_plus_centroidalifold_bpp_mat_file_path + " " + mafft_xinsi_output_file_path
    mafft_xinsi_plus_centroidalifold_params_4_bpp_mat.insert(0, mafft_xinsi_plus_centroidalifold_command)
    ref_sa_plus_centroidalifold_command = "centroid_alifold -e McCaskill -w 0 -e Alifold -w 1 --posteriors 0 --posteriors-output " + ref_sa_plus_centroidalifold_bpp_mat_file_path + " " + ref_sa_file_path
    ref_sa_plus_centroidalifold_params_4_bpp_mat.insert(0, ref_sa_plus_centroidalifold_command)
  pool = multiprocessing.Pool(num_of_threads)
  pool.map(utils.run_command, mafft_ginsi_plus_centroidalifold_params_4_bpp_mat)
  begin = time.time()
  pool.map(utils.run_command, mafft_xinsi_plus_centroidalifold_params_4_bpp_mat)
  phyloalifold_elapsed_time += time.time() - begin
  pool.map(utils.run_command, ref_sa_plus_centroidalifold_params_4_bpp_mat)
  sub_thread_num = 4
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
    mafft_ginsi_plus_phyloalifold_output_dir_path = os.path.join(mafft_ginsi_plus_phyloalifold_dir_path, "csss_of_" + rna_family_name)
    mafft_xinsi_plus_phyloalifold_output_dir_path = os.path.join(mafft_xinsi_plus_phyloalifold_dir_path, "csss_of_" + rna_family_name)
    ref_sa_plus_phyloalifold_output_dir_path = os.path.join(ref_sa_plus_phyloalifold_dir_path, "csss_of_" + rna_family_name)
    mafft_ginsi_plus_centroidalifold_output_dir_path = os.path.join(mafft_ginsi_plus_centroidalifold_dir_path, "csss_of_" + rna_family_name)
    mafft_xinsi_plus_centroidalifold_output_dir_path = os.path.join(mafft_xinsi_plus_centroidalifold_dir_path, "csss_of_" + rna_family_name)
    ref_sa_plus_centroidalifold_output_dir_path = os.path.join(ref_sa_plus_centroidalifold_dir_path, "csss_of_" + rna_family_name)
    if not os.path.isdir(mafft_ginsi_plus_phyloalifold_output_dir_path):
      os.mkdir(mafft_ginsi_plus_phyloalifold_output_dir_path)
    if not os.path.isdir(mafft_xinsi_plus_phyloalifold_output_dir_path):
      os.mkdir(mafft_xinsi_plus_phyloalifold_output_dir_path)
    if not os.path.isdir(ref_sa_plus_phyloalifold_output_dir_path):
      os.mkdir(ref_sa_plus_phyloalifold_output_dir_path)
    if not os.path.isdir(mafft_ginsi_plus_centroidalifold_output_dir_path):
      os.mkdir(mafft_ginsi_plus_centroidalifold_output_dir_path)
    if not os.path.isdir(mafft_xinsi_plus_centroidalifold_output_dir_path):
      os.mkdir(mafft_xinsi_plus_centroidalifold_output_dir_path)
    if not os.path.isdir(ref_sa_plus_centroidalifold_output_dir_path):
      os.mkdir(ref_sa_plus_centroidalifold_output_dir_path)
    phyloalifold_command = "phyloalifold -t " + str(sub_thread_num) + " -i " + rna_seq_file_path + " -a " + mafft_ginsi_output_file_path + " -c " + mafft_ginsi_plus_centroidalifold_bpp_mat_file_path + " -o " + mafft_ginsi_plus_phyloalifold_output_dir_path
    mafft_ginsi_plus_phyloalifold_params.insert(0, phyloalifold_command)
    mafft_ginsi_plus_centroidalifold_bpp_mat_file_path = os.path.join(mafft_ginsi_plus_centroidalifold_bpp_mat_dir_path, rna_family_name + ".dat")
    phyloalifold_command = "phyloalifold -t " + str(sub_thread_num) + " -i " + rna_seq_file_path + " -a " + mafft_xinsi_output_file_path + " -c " + mafft_xinsi_plus_centroidalifold_bpp_mat_file_path + " -o " + mafft_xinsi_plus_phyloalifold_output_dir_path
    mafft_xinsi_plus_phyloalifold_params.insert(0, phyloalifold_command)
    phyloalifold_command = "phyloalifold -t " + str(sub_thread_num) + " -i " + rna_seq_file_path + " -a " + ref_sa_file_path + " -c " + ref_sa_plus_centroidalifold_bpp_mat_file_path + " -o " + ref_sa_plus_phyloalifold_output_dir_path
    ref_sa_plus_phyloalifold_params.insert(0, phyloalifold_command)
    for gamma in gammas:
      gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
      output_file = "gamma=" + gamma_str + ".sth"
      mafft_ginsi_plus_centroidalifold_output_file_path = os.path.join(mafft_ginsi_plus_centroidalifold_output_dir_path, output_file)
      mafft_xinsi_plus_centroidalifold_output_file_path = os.path.join(mafft_xinsi_plus_centroidalifold_output_dir_path, output_file)
      ref_sa_plus_centroidalifold_output_file_path = os.path.join(ref_sa_plus_centroidalifold_output_dir_path, output_file)
      mafft_ginsi_plus_centroidalifold_params.insert(0, (mafft_ginsi_output_file_path, mafft_ginsi_plus_centroidalifold_output_file_path, gamma_str))
      mafft_xinsi_plus_centroidalifold_params.insert(0, (mafft_xinsi_output_file_path, mafft_xinsi_plus_centroidalifold_output_file_path, gamma_str))
      ref_sa_plus_centroidalifold_params.insert(0, (ref_sa_file_path, ref_sa_plus_centroidalifold_output_file_path, gamma_str))
      if gamma == 1:
        centroidalifold_params_4_elapsed_time.insert(0, (mafft_xinsi_output_file_path, mafft_xinsi_plus_centroidalifold_output_file_path, gamma_str))
  pool = multiprocessing.Pool(int(num_of_threads / sub_thread_num))
  pool.map(utils.run_command, mafft_ginsi_plus_phyloalifold_params)
  begin = time.time()
  pool.map(utils.run_command, mafft_xinsi_plus_phyloalifold_params)
  phyloalifold_elapsed_time += time.time() - begin
  pool.map(utils.run_command, ref_sa_plus_phyloalifold_params)
  pool = multiprocessing.Pool(num_of_threads)
  pool.map(run_centroidalifold, mafft_ginsi_plus_centroidalifold_params)
  pool.map(run_centroidalifold, mafft_xinsi_plus_centroidalifold_params)
  pool.map(run_centroidalifold, ref_sa_plus_centroidalifold_params)
  begin = time.time()
  pool.map(run_centroidalifold, centroidalifold_params_4_elapsed_time)
  centroidalifold_elapsed_time = time.time() - begin
  print("The elapsed time of the PhyloAliFold program for a test set = %f [s]." % phyloalifold_elapsed_time)
  print("The elapsed time of the CentroidAlifold program for a test set = %f [s]." % centroidalifold_elapsed_time)

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

if __name__ == "__main__":
  main()
