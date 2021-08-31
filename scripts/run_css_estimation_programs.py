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
  gammas = [2. ** i for i in range(-4, 11)]
  mafft_params = []
  probcons_params = []
  clustalw_params = []
  mafft_xinsi_params = []
  mafft_plus_consalifold_params = []
  probcons_plus_consalifold_params = []
  clustalw_plus_consalifold_params = []
  mafft_xinsi_plus_consalifold_params = []
  ref_sa_plus_consalifold_params = []
  posterior_probcons_plus_consalifold_params = []
  posterior_clustalw_plus_consalifold_params = []
  posterior_mafft_plus_consalifold_params = []
  posterior_mafft_xinsi_plus_consalifold_params = []
  posterior_ref_sa_plus_consalifold_params = []
  contra_mafft_plus_consalifold_params = []
  contra_probcons_plus_consalifold_params = []
  contra_clustalw_plus_consalifold_params = []
  contra_mafft_xinsi_plus_consalifold_params = []
  contra_ref_sa_plus_consalifold_params = []
  consalifold_params_4_elapsed_time = []
  mafft_plus_centroidalifold_params = []
  probcons_plus_centroidalifold_params = []
  clustalw_plus_centroidalifold_params = []
  mafft_xinsi_plus_centroidalifold_params = []
  ref_sa_plus_centroidalifold_params = []
  centroidalifold_params_4_elapsed_time = []
  mafft_plus_rnaalifold_params = []
  probcons_plus_rnaalifold_params = []
  clustalw_plus_rnaalifold_params = []
  mafft_xinsi_plus_rnaalifold_params = []
  ref_sa_plus_rnaalifold_params = []
  mafft_plus_petfold_params = []
  probcons_plus_petfold_params = []
  clustalw_plus_petfold_params = []
  mafft_xinsi_plus_petfold_params = []
  ref_sa_plus_petfold_params = []
  petfold_params_4_elapsed_time = []
  rna_seq_dir_path = asset_dir_path + "/compiled_rna_fams_test"
  mafft_dir_path = asset_dir_path + "/mafft"
  probcons_dir_path = asset_dir_path + "/probcons"
  clustalw_dir_path = asset_dir_path + "/clustalw"
  mafft_xinsi_dir_path = asset_dir_path + "/mafft_xinsi"
  ref_sa_dir_path = asset_dir_path + "/ref_sas_test"
  mafft_plus_consalifold_dir_path = asset_dir_path + "/mafft_plus_consalifold"
  probcons_plus_consalifold_dir_path = asset_dir_path + "/probcons_plus_consalifold"
  clustalw_plus_consalifold_dir_path = asset_dir_path + "/clustalw_plus_consalifold"
  mafft_xinsi_plus_consalifold_dir_path = asset_dir_path + "/mafft_xinsi_plus_consalifold"
  ref_sa_plus_consalifold_dir_path = asset_dir_path + "/ref_sa_plus_consalifold"
  contra_mafft_plus_consalifold_dir_path = asset_dir_path + "/contra_mafft_plus_consalifold"
  contra_probcons_plus_consalifold_dir_path = asset_dir_path + "/contra_probcons_plus_consalifold"
  contra_clustalw_plus_consalifold_dir_path = asset_dir_path + "/contra_clustalw_plus_consalifold"
  contra_mafft_xinsi_plus_consalifold_dir_path = asset_dir_path + "/contra_mafft_xinsi_plus_consalifold"
  contra_ref_sa_plus_consalifold_dir_path = asset_dir_path + "/contra_ref_sa_plus_consalifold"
  posterior_probcons_plus_consalifold_dir_path = asset_dir_path + "/posterior_probcons_plus_consalifold"
  posterior_clustalw_plus_consalifold_dir_path = asset_dir_path + "/posterior_clustalw_plus_consalifold"
  posterior_mafft_plus_consalifold_dir_path = asset_dir_path + "/posterior_mafft_plus_consalifold"
  posterior_mafft_xinsi_plus_consalifold_dir_path = asset_dir_path + "/posterior_mafft_xinsi_plus_consalifold"
  posterior_ref_sa_plus_consalifold_dir_path = asset_dir_path + "/posterior_ref_sa_plus_consalifold"
  mafft_plus_centroidalifold_dir_path = asset_dir_path + "/mafft_plus_centroidalifold"
  probcons_plus_centroidalifold_dir_path = asset_dir_path + "/probcons_plus_centroidalifold"
  clustalw_plus_centroidalifold_dir_path = asset_dir_path + "/clustalw_plus_centroidalifold"
  mafft_xinsi_plus_centroidalifold_dir_path = asset_dir_path + "/mafft_xinsi_plus_centroidalifold"
  ref_sa_plus_centroidalifold_dir_path = asset_dir_path + "/ref_sa_plus_centroidalifold"
  mafft_plus_rnaalifold_dir_path = asset_dir_path + "/mafft_plus_rnaalifold"
  probcons_plus_rnaalifold_dir_path = asset_dir_path + "/probcons_plus_rnaalifold"
  clustalw_plus_rnaalifold_dir_path = asset_dir_path + "/clustalw_plus_rnaalifold"
  mafft_xinsi_plus_rnaalifold_dir_path = asset_dir_path + "/mafft_xinsi_plus_rnaalifold"
  ref_sa_plus_rnaalifold_dir_path = asset_dir_path + "/ref_sa_plus_rnaalifold"
  mafft_plus_petfold_dir_path = asset_dir_path + "/mafft_plus_petfold"
  probcons_plus_petfold_dir_path = asset_dir_path + "/probcons_plus_petfold"
  clustalw_plus_petfold_dir_path = asset_dir_path + "/clustalw_plus_petfold"
  mafft_xinsi_plus_petfold_dir_path = asset_dir_path + "/mafft_xinsi_plus_petfold"
  ref_sa_plus_petfold_dir_path = asset_dir_path + "/ref_sa_plus_petfold"
  if not os.path.isdir(mafft_dir_path):
    os.mkdir(mafft_dir_path)
  if not os.path.isdir(probcons_dir_path):
    os.mkdir(probcons_dir_path)
  if not os.path.isdir(clustalw_dir_path):
    os.mkdir(clustalw_dir_path)
  if not os.path.isdir(mafft_xinsi_dir_path):
    os.mkdir(mafft_xinsi_dir_path)
  if not os.path.isdir(mafft_plus_consalifold_dir_path):
    os.mkdir(mafft_plus_consalifold_dir_path)
  if not os.path.isdir(probcons_plus_consalifold_dir_path):
    os.mkdir(probcons_plus_consalifold_dir_path)
  if not os.path.isdir(clustalw_plus_consalifold_dir_path):
    os.mkdir(clustalw_plus_consalifold_dir_path)
  if not os.path.isdir(mafft_xinsi_plus_consalifold_dir_path):
    os.mkdir(mafft_xinsi_plus_consalifold_dir_path)
  if not os.path.isdir(ref_sa_plus_consalifold_dir_path):
    os.mkdir(ref_sa_plus_consalifold_dir_path)
  if not os.path.isdir(contra_mafft_plus_consalifold_dir_path):
    os.mkdir(contra_mafft_plus_consalifold_dir_path)
  if not os.path.isdir(contra_probcons_plus_consalifold_dir_path):
    os.mkdir(contra_probcons_plus_consalifold_dir_path)
  if not os.path.isdir(contra_clustalw_plus_consalifold_dir_path):
    os.mkdir(contra_clustalw_plus_consalifold_dir_path)
  if not os.path.isdir(contra_mafft_xinsi_plus_consalifold_dir_path):
    os.mkdir(contra_mafft_xinsi_plus_consalifold_dir_path)
  if not os.path.isdir(contra_ref_sa_plus_consalifold_dir_path):
    os.mkdir(contra_ref_sa_plus_consalifold_dir_path)
  if not os.path.isdir(posterior_probcons_plus_consalifold_dir_path):
    os.mkdir(posterior_probcons_plus_consalifold_dir_path)
  if not os.path.isdir(posterior_clustalw_plus_consalifold_dir_path):
    os.mkdir(posterior_clustalw_plus_consalifold_dir_path)
  if not os.path.isdir(posterior_mafft_plus_consalifold_dir_path):
    os.mkdir(posterior_mafft_plus_consalifold_dir_path)
  if not os.path.isdir(posterior_mafft_xinsi_plus_consalifold_dir_path):
    os.mkdir(posterior_mafft_xinsi_plus_consalifold_dir_path)
  if not os.path.isdir(posterior_ref_sa_plus_consalifold_dir_path):
    os.mkdir(posterior_ref_sa_plus_consalifold_dir_path)
  if not os.path.isdir(mafft_plus_centroidalifold_dir_path):
    os.mkdir(mafft_plus_centroidalifold_dir_path)
  if not os.path.isdir(probcons_plus_centroidalifold_dir_path):
    os.mkdir(probcons_plus_centroidalifold_dir_path)
  if not os.path.isdir(clustalw_plus_centroidalifold_dir_path):
    os.mkdir(clustalw_plus_centroidalifold_dir_path)
  if not os.path.isdir(mafft_xinsi_plus_centroidalifold_dir_path):
    os.mkdir(mafft_xinsi_plus_centroidalifold_dir_path)
  if not os.path.isdir(ref_sa_plus_centroidalifold_dir_path):
    os.mkdir(ref_sa_plus_centroidalifold_dir_path)
  if not os.path.isdir(mafft_plus_rnaalifold_dir_path):
    os.mkdir(mafft_plus_rnaalifold_dir_path)
  if not os.path.isdir(probcons_plus_rnaalifold_dir_path):
    os.mkdir(probcons_plus_rnaalifold_dir_path)
  if not os.path.isdir(clustalw_plus_rnaalifold_dir_path):
    os.mkdir(clustalw_plus_rnaalifold_dir_path)
  if not os.path.isdir(mafft_xinsi_plus_rnaalifold_dir_path):
    os.mkdir(mafft_xinsi_plus_rnaalifold_dir_path)
  if not os.path.isdir(ref_sa_plus_rnaalifold_dir_path):
    os.mkdir(ref_sa_plus_rnaalifold_dir_path)
  if not os.path.isdir(mafft_plus_petfold_dir_path):
    os.mkdir(mafft_plus_petfold_dir_path)
  if not os.path.isdir(probcons_plus_petfold_dir_path):
    os.mkdir(probcons_plus_petfold_dir_path)
  if not os.path.isdir(clustalw_plus_petfold_dir_path):
    os.mkdir(clustalw_plus_petfold_dir_path)
  if not os.path.isdir(mafft_xinsi_plus_petfold_dir_path):
    os.mkdir(mafft_xinsi_plus_petfold_dir_path)
  if not os.path.isdir(ref_sa_plus_petfold_dir_path):
    os.mkdir(ref_sa_plus_petfold_dir_path)
  for rna_seq_file in os.listdir(rna_seq_dir_path):
    if not rna_seq_file.endswith(".fa"):
      continue
    rna_seq_file_path = os.path.join(rna_seq_dir_path, rna_seq_file)
    (rna_family_name, extension) = os.path.splitext(rna_seq_file)
    mafft_output_file_path = os.path.join(mafft_dir_path, rna_family_name + ".aln")
    probcons_output_file_path = os.path.join(probcons_dir_path, rna_family_name + ".aln")
    clustalw_output_file_path = os.path.join(clustalw_dir_path, rna_family_name + ".aln")
    mafft_xinsi_output_file_path = os.path.join(mafft_xinsi_dir_path, rna_family_name + ".aln")
    ref_sa_file_path = os.path.join(ref_sa_dir_path, rna_family_name + ".aln")
    mafft_params.insert(0, (rna_seq_file_path, mafft_output_file_path))
    probcons_params.insert(0, (rna_seq_file_path, probcons_output_file_path))
    clustalw_params.insert(0, (rna_seq_file_path, clustalw_output_file_path))
    mafft_xinsi_params.insert(0, (rna_seq_file_path, mafft_xinsi_output_file_path))
  pool = multiprocessing.Pool(num_of_threads)
  pool.map(run_mafft, mafft_params)
  pool.map(run_probcons, probcons_params)
  pool.map(run_clustalw, clustalw_params)
  pool.map(run_mafft_xinsi, mafft_xinsi_params)
  sub_thread_num = 4
  for rna_seq_file in os.listdir(rna_seq_dir_path):
    if not rna_seq_file.endswith(".fa"):
      continue
    rna_seq_file_path = os.path.join(rna_seq_dir_path, rna_seq_file)
    (rna_family_name, extension) = os.path.splitext(rna_seq_file)
    mafft_output_file_path = os.path.join(mafft_dir_path, rna_family_name + ".aln")
    probcons_output_file_path = os.path.join(probcons_dir_path, rna_family_name + ".aln")
    clustalw_output_file_path = os.path.join(clustalw_dir_path, rna_family_name + ".aln")
    mafft_xinsi_output_file_path = os.path.join(mafft_xinsi_dir_path, rna_family_name + ".aln")
    ref_sa_file_path = os.path.join(ref_sa_dir_path, rna_family_name + ".aln")
    mafft_plus_consalifold_output_dir_path = os.path.join(mafft_plus_consalifold_dir_path, rna_family_name)
    probcons_plus_consalifold_output_dir_path = os.path.join(probcons_plus_consalifold_dir_path, rna_family_name)
    clustalw_plus_consalifold_output_dir_path = os.path.join(clustalw_plus_consalifold_dir_path, rna_family_name)
    mafft_xinsi_plus_consalifold_output_dir_path = os.path.join(mafft_xinsi_plus_consalifold_dir_path, rna_family_name)
    ref_sa_plus_consalifold_output_dir_path = os.path.join(ref_sa_plus_consalifold_dir_path, rna_family_name)
    contra_mafft_plus_consalifold_output_dir_path = os.path.join(contra_mafft_plus_consalifold_dir_path, rna_family_name)
    contra_probcons_plus_consalifold_output_dir_path = os.path.join(contra_probcons_plus_consalifold_dir_path, rna_family_name)
    contra_clustalw_plus_consalifold_output_dir_path = os.path.join(contra_clustalw_plus_consalifold_dir_path, rna_family_name)
    contra_mafft_xinsi_plus_consalifold_output_dir_path = os.path.join(contra_mafft_xinsi_plus_consalifold_dir_path, rna_family_name)
    contra_ref_sa_plus_consalifold_output_dir_path = os.path.join(contra_ref_sa_plus_consalifold_dir_path, rna_family_name)
    posterior_probcons_plus_consalifold_output_dir_path = os.path.join(posterior_probcons_plus_consalifold_dir_path, rna_family_name)
    posterior_clustalw_plus_consalifold_output_dir_path = os.path.join(posterior_clustalw_plus_consalifold_dir_path, rna_family_name)
    posterior_mafft_plus_consalifold_output_dir_path = os.path.join(posterior_mafft_plus_consalifold_dir_path, rna_family_name)
    posterior_mafft_xinsi_plus_consalifold_output_dir_path = os.path.join(posterior_mafft_xinsi_plus_consalifold_dir_path, rna_family_name)
    posterior_ref_sa_plus_consalifold_output_dir_path = os.path.join(posterior_ref_sa_plus_consalifold_dir_path, rna_family_name)
    mafft_plus_centroidalifold_output_dir_path = os.path.join(mafft_plus_centroidalifold_dir_path, rna_family_name)
    probcons_plus_centroidalifold_output_dir_path = os.path.join(probcons_plus_centroidalifold_dir_path, rna_family_name)
    clustalw_plus_centroidalifold_output_dir_path = os.path.join(clustalw_plus_centroidalifold_dir_path, rna_family_name)
    mafft_xinsi_plus_centroidalifold_output_dir_path = os.path.join(mafft_xinsi_plus_centroidalifold_dir_path, rna_family_name)
    ref_sa_plus_centroidalifold_output_dir_path = os.path.join(ref_sa_plus_centroidalifold_dir_path, rna_family_name)
    mafft_plus_petfold_output_dir_path = os.path.join(mafft_plus_petfold_dir_path, rna_family_name)
    probcons_plus_petfold_output_dir_path = os.path.join(probcons_plus_petfold_dir_path, rna_family_name)
    clustalw_plus_petfold_output_dir_path = os.path.join(clustalw_plus_petfold_dir_path, rna_family_name)
    mafft_xinsi_plus_petfold_output_dir_path = os.path.join(mafft_xinsi_plus_petfold_dir_path, rna_family_name)
    ref_sa_plus_petfold_output_dir_path = os.path.join(ref_sa_plus_petfold_dir_path, rna_family_name)
    if not os.path.isdir(mafft_plus_consalifold_output_dir_path):
      os.mkdir(mafft_plus_consalifold_output_dir_path)
    if not os.path.isdir(probcons_plus_consalifold_output_dir_path):
      os.mkdir(probcons_plus_consalifold_output_dir_path)
    if not os.path.isdir(posterior_probcons_plus_consalifold_output_dir_path):
      os.mkdir(posterior_probcons_plus_consalifold_output_dir_path)
    if not os.path.isdir(posterior_clustalw_plus_consalifold_output_dir_path):
      os.mkdir(posterior_clustalw_plus_consalifold_output_dir_path)
    if not os.path.isdir(posterior_mafft_plus_consalifold_output_dir_path):
      os.mkdir(posterior_mafft_plus_consalifold_output_dir_path)
    if not os.path.isdir(posterior_mafft_xinsi_plus_consalifold_output_dir_path):
      os.mkdir(posterior_mafft_xinsi_plus_consalifold_output_dir_path)
    if not os.path.isdir(posterior_ref_sa_plus_consalifold_output_dir_path):
      os.mkdir(posterior_ref_sa_plus_consalifold_output_dir_path)
    if not os.path.isdir(clustalw_plus_consalifold_output_dir_path):
      os.mkdir(clustalw_plus_consalifold_output_dir_path)
    if not os.path.isdir(mafft_xinsi_plus_consalifold_output_dir_path):
      os.mkdir(mafft_xinsi_plus_consalifold_output_dir_path)
    if not os.path.isdir(ref_sa_plus_consalifold_output_dir_path):
      os.mkdir(ref_sa_plus_consalifold_output_dir_path)
    if not os.path.isdir(mafft_plus_centroidalifold_output_dir_path):
      os.mkdir(mafft_plus_centroidalifold_output_dir_path)
    if not os.path.isdir(probcons_plus_centroidalifold_output_dir_path):
      os.mkdir(probcons_plus_centroidalifold_output_dir_path)
    if not os.path.isdir(clustalw_plus_centroidalifold_output_dir_path):
      os.mkdir(clustalw_plus_centroidalifold_output_dir_path)
    if not os.path.isdir(mafft_xinsi_plus_centroidalifold_output_dir_path):
      os.mkdir(mafft_xinsi_plus_centroidalifold_output_dir_path)
    if not os.path.isdir(ref_sa_plus_centroidalifold_output_dir_path):
      os.mkdir(ref_sa_plus_centroidalifold_output_dir_path)
    if not os.path.isdir(mafft_plus_petfold_output_dir_path):
      os.mkdir(mafft_plus_petfold_output_dir_path)
    if not os.path.isdir(probcons_plus_petfold_output_dir_path):
      os.mkdir(probcons_plus_petfold_output_dir_path)
    if not os.path.isdir(clustalw_plus_petfold_output_dir_path):
      os.mkdir(clustalw_plus_petfold_output_dir_path)
    if not os.path.isdir(mafft_xinsi_plus_petfold_output_dir_path):
      os.mkdir(mafft_xinsi_plus_petfold_output_dir_path)
    if not os.path.isdir(ref_sa_plus_petfold_output_dir_path):
      os.mkdir(ref_sa_plus_petfold_output_dir_path)
    mafft_plus_consalifold_command = "consalifold -t " + str(sub_thread_num) + " -i " + mafft_output_file_path + " -o " + mafft_plus_consalifold_output_dir_path
    probcons_plus_consalifold_command = "consalifold -t " + str(sub_thread_num) + " -i " + probcons_output_file_path + " -o " + probcons_plus_consalifold_output_dir_path
    clustalw_plus_consalifold_command = "consalifold -t " + str(sub_thread_num) + " -i " + clustalw_output_file_path + " -o " + clustalw_plus_consalifold_output_dir_path
    mafft_xinsi_plus_consalifold_command = "consalifold -t " + str(sub_thread_num) + " -i " + mafft_xinsi_output_file_path + " -o " + mafft_xinsi_plus_consalifold_output_dir_path
    ref_sa_plus_consalifold_command = "consalifold -t " + str(sub_thread_num) + " -i " + ref_sa_file_path + " -o " + ref_sa_plus_consalifold_output_dir_path
    mafft_plus_consalifold_params.insert(0, mafft_plus_consalifold_command)
    probcons_plus_consalifold_params.insert(0, probcons_plus_consalifold_command)
    clustalw_plus_consalifold_params.insert(0, clustalw_plus_consalifold_command)
    mafft_xinsi_plus_consalifold_params.insert(0, mafft_xinsi_plus_consalifold_command)
    ref_sa_plus_consalifold_params.insert(0, ref_sa_plus_consalifold_command)
    mafft_xinsi_plus_consalifold_command = "consalifold -g 1 -t " + str(sub_thread_num) + " -i " + mafft_xinsi_output_file_path + " -o " + mafft_xinsi_plus_consalifold_output_dir_path
    consalifold_params_4_elapsed_time.insert(0, mafft_xinsi_plus_consalifold_command)
    contra_mafft_plus_consalifold_command = "consalifold -m contra -t " + str(sub_thread_num) + " -i " + mafft_output_file_path + " -o " + contra_mafft_plus_consalifold_output_dir_path
    contra_probcons_plus_consalifold_command = "consalifold -m contra -t " + str(sub_thread_num) + " -i " + probcons_output_file_path + " -o " + contra_probcons_plus_consalifold_output_dir_path
    contra_clustalw_plus_consalifold_command = "consalifold -m contra -t " + str(sub_thread_num) + " -i " + clustalw_output_file_path + " -o " + contra_clustalw_plus_consalifold_output_dir_path
    contra_mafft_xinsi_plus_consalifold_command = "consalifold -m contra -t " + str(sub_thread_num) + " -i " + mafft_xinsi_output_file_path + " -o " + contra_mafft_xinsi_plus_consalifold_output_dir_path
    contra_ref_sa_plus_consalifold_command = "consalifold -m contra -t " + str(sub_thread_num) + " -i " + ref_sa_file_path + " -o " + contra_ref_sa_plus_consalifold_output_dir_path
    contra_mafft_plus_consalifold_params.insert(0, contra_mafft_plus_consalifold_command)
    contra_probcons_plus_consalifold_params.insert(0, contra_probcons_plus_consalifold_command)
    contra_clustalw_plus_consalifold_params.insert(0, contra_clustalw_plus_consalifold_command)
    contra_mafft_xinsi_plus_consalifold_params.insert(0, contra_mafft_xinsi_plus_consalifold_command)
    contra_ref_sa_plus_consalifold_params.insert(0, contra_ref_sa_plus_consalifold_command)
    posterior_probcons_plus_consalifold_command = "consalifold -m posterior -t " + str(sub_thread_num) + " -i " + probcons_output_file_path + " -o " + posterior_probcons_plus_consalifold_output_dir_path
    posterior_probcons_plus_consalifold_params.insert(0, posterior_probcons_plus_consalifold_command)
    posterior_clustalw_plus_consalifold_command = "consalifold -m posterior -t " + str(sub_thread_num) + " -i " + clustalw_output_file_path + " -o " + posterior_clustalw_plus_consalifold_output_dir_path
    posterior_clustalw_plus_consalifold_params.insert(0, posterior_clustalw_plus_consalifold_command)
    posterior_mafft_plus_consalifold_command = "consalifold -m posterior -t " + str(sub_thread_num) + " -i " + mafft_output_file_path + " -o " + posterior_mafft_plus_consalifold_output_dir_path
    posterior_mafft_plus_consalifold_params.insert(0, posterior_mafft_plus_consalifold_command)
    posterior_mafft_xinsi_plus_consalifold_command = "consalifold -m posterior -t " + str(sub_thread_num) + " -i " + mafft_xinsi_output_file_path + " -o " + posterior_mafft_xinsi_plus_consalifold_output_dir_path
    posterior_mafft_xinsi_plus_consalifold_params.insert(0, posterior_mafft_xinsi_plus_consalifold_command)
    posterior_ref_sa_plus_consalifold_command = "consalifold -m posterior -t " + str(sub_thread_num) + " -i " + ref_sa_file_path + " -o " + posterior_ref_sa_plus_consalifold_output_dir_path
    posterior_ref_sa_plus_consalifold_params.insert(0, posterior_ref_sa_plus_consalifold_command)
    output_file = rna_family_name + ".sth"
    mafft_plus_rnaalifold_output_file_path = os.path.join(mafft_plus_rnaalifold_dir_path, output_file)
    probcons_plus_rnaalifold_output_file_path = os.path.join(probcons_plus_rnaalifold_dir_path, output_file)
    clustalw_plus_rnaalifold_output_file_path = os.path.join(clustalw_plus_rnaalifold_dir_path, output_file)
    mafft_xinsi_plus_rnaalifold_output_file_path = os.path.join(mafft_xinsi_plus_rnaalifold_dir_path, output_file)
    ref_sa_plus_rnaalifold_output_file_path = os.path.join(ref_sa_plus_rnaalifold_dir_path, output_file)
    mafft_plus_rnaalifold_params.insert(0, (mafft_output_file_path, mafft_plus_rnaalifold_output_file_path))
    probcons_plus_rnaalifold_params.insert(0, (probcons_output_file_path, probcons_plus_rnaalifold_output_file_path))
    clustalw_plus_rnaalifold_params.insert(0, (clustalw_output_file_path, clustalw_plus_rnaalifold_output_file_path))
    mafft_xinsi_plus_rnaalifold_params.insert(0, (mafft_xinsi_output_file_path, mafft_xinsi_plus_rnaalifold_output_file_path))
    ref_sa_plus_rnaalifold_params.insert(0, (ref_sa_file_path, ref_sa_plus_rnaalifold_output_file_path))
    mafft_output_file_path_2 = os.path.join(mafft_dir_path, rna_family_name + ".fa")
    probcons_output_file_path_2 = os.path.join(probcons_dir_path, rna_family_name + ".fa")
    clustalw_output_file_path_2 = os.path.join(clustalw_dir_path, rna_family_name + ".fa")
    mafft_xinsi_output_file_path_2 = os.path.join(mafft_xinsi_dir_path, rna_family_name + ".fa")
    ref_sa_file_path_2 = os.path.join(ref_sa_dir_path, rna_family_name + ".fa")
    for gamma in gammas:
      gamma_str = str(gamma) if gamma < 1 else str(int(gamma))
      output_file = "gamma=" + gamma_str + ".sth"
      petfold_gamma_str = str(1 / gamma) if gamma > 1 else str(int(1 / gamma))
      mafft_plus_centroidalifold_output_file_path = os.path.join(mafft_plus_centroidalifold_output_dir_path, output_file)
      probcons_plus_centroidalifold_output_file_path = os.path.join(probcons_plus_centroidalifold_output_dir_path, output_file)
      clustalw_plus_centroidalifold_output_file_path = os.path.join(clustalw_plus_centroidalifold_output_dir_path, output_file)
      mafft_xinsi_plus_centroidalifold_output_file_path = os.path.join(mafft_xinsi_plus_centroidalifold_output_dir_path, output_file)
      ref_sa_plus_centroidalifold_output_file_path = os.path.join(ref_sa_plus_centroidalifold_output_dir_path, output_file)
      mafft_plus_centroidalifold_params.insert(0, (mafft_output_file_path, mafft_plus_centroidalifold_output_file_path, gamma_str))
      probcons_plus_centroidalifold_params.insert(0, (probcons_output_file_path, probcons_plus_centroidalifold_output_file_path, gamma_str))
      clustalw_plus_centroidalifold_params.insert(0, (clustalw_output_file_path, clustalw_plus_centroidalifold_output_file_path, gamma_str))
      mafft_xinsi_plus_centroidalifold_params.insert(0, (mafft_xinsi_output_file_path, mafft_xinsi_plus_centroidalifold_output_file_path, gamma_str))
      ref_sa_plus_centroidalifold_params.insert(0, (ref_sa_file_path, ref_sa_plus_centroidalifold_output_file_path, gamma_str))
      if gamma == 1:
        centroidalifold_params_4_elapsed_time.insert(0, (mafft_xinsi_output_file_path, mafft_xinsi_plus_centroidalifold_output_file_path, gamma_str))
      mafft_plus_petfold_output_file_path = os.path.join(mafft_plus_petfold_output_dir_path, output_file)
      probcons_plus_petfold_output_file_path = os.path.join(probcons_plus_petfold_output_dir_path, output_file)
      clustalw_plus_petfold_output_file_path = os.path.join(clustalw_plus_petfold_output_dir_path, output_file)
      mafft_xinsi_plus_petfold_output_file_path = os.path.join(mafft_xinsi_plus_petfold_output_dir_path, output_file)
      ref_sa_plus_petfold_output_file_path = os.path.join(ref_sa_plus_petfold_output_dir_path, output_file)
      mafft_plus_petfold_params.insert(0, (mafft_output_file_path_2, mafft_plus_petfold_output_file_path, petfold_gamma_str))
      probcons_plus_petfold_params.insert(0, (probcons_output_file_path_2, probcons_plus_petfold_output_file_path, petfold_gamma_str))
      clustalw_plus_petfold_params.insert(0, (clustalw_output_file_path_2, clustalw_plus_petfold_output_file_path, petfold_gamma_str))
      mafft_xinsi_plus_petfold_params.insert(0, (mafft_xinsi_output_file_path_2, mafft_xinsi_plus_petfold_output_file_path, petfold_gamma_str))
      ref_sa_plus_petfold_params.insert(0, (ref_sa_file_path_2, ref_sa_plus_petfold_output_file_path, petfold_gamma_str))
      if gamma == 1:
        petfold_params_4_elapsed_time.insert(0, (mafft_xinsi_output_file_path_2, mafft_xinsi_plus_petfold_output_file_path, petfold_gamma_str))
  # ConsAliFold's execution.
  pool = multiprocessing.Pool(int(num_of_threads / sub_thread_num))
  pool.map(utils.run_command, mafft_plus_consalifold_params)
  pool.map(utils.run_command, probcons_plus_consalifold_params)
  pool.map(utils.run_command, clustalw_plus_consalifold_params)
  pool.map(utils.run_command, mafft_xinsi_plus_consalifold_params)
  pool.map(utils.run_command, ref_sa_plus_consalifold_params)
  pool.map(utils.run_command, contra_mafft_plus_consalifold_params)
  pool.map(utils.run_command, contra_probcons_plus_consalifold_params)
  pool.map(utils.run_command, contra_clustalw_plus_consalifold_params)
  pool.map(utils.run_command, contra_mafft_xinsi_plus_consalifold_params)
  pool.map(utils.run_command, contra_ref_sa_plus_consalifold_params)
  pool.map(utils.run_command, posterior_probcons_plus_consalifold_params)
  pool.map(utils.run_command, posterior_clustalw_plus_consalifold_params)
  pool.map(utils.run_command, posterior_mafft_plus_consalifold_params)
  pool.map(utils.run_command, posterior_mafft_xinsi_plus_consalifold_params)
  pool.map(utils.run_command, posterior_ref_sa_plus_consalifold_params)
  begin = time.time()
  pool.map(utils.run_command, consalifold_params_4_elapsed_time)
  consalifold_elapsed_time = time.time() - begin
  # CentroidAlifold's execution.
  pool = multiprocessing.Pool(num_of_threads)
  pool.map(run_centroidalifold, mafft_plus_centroidalifold_params)
  pool.map(run_centroidalifold, probcons_plus_centroidalifold_params)
  pool.map(run_centroidalifold, clustalw_plus_centroidalifold_params)
  pool.map(run_centroidalifold, mafft_xinsi_plus_centroidalifold_params)
  pool.map(run_centroidalifold, ref_sa_plus_centroidalifold_params)
  begin = time.time()
  pool.map(run_centroidalifold, centroidalifold_params_4_elapsed_time)
  centroidalifold_elapsed_time = time.time() - begin
  # PetFold's execution.
  pool.map(run_petfold, mafft_plus_petfold_params)
  pool.map(run_petfold, probcons_plus_petfold_params)
  pool.map(run_petfold, clustalw_plus_petfold_params)
  pool.map(run_petfold, mafft_xinsi_plus_petfold_params)
  pool.map(run_petfold, ref_sa_plus_petfold_params)
  begin = time.time()
  pool.map(run_petfold, petfold_params_4_elapsed_time)
  petfold_elapsed_time = time.time() - begin
  # RNAalifold's execution.
  pool.map(run_rnaalifold, mafft_plus_rnaalifold_params)
  pool.map(run_rnaalifold, probcons_plus_rnaalifold_params)
  pool.map(run_rnaalifold, clustalw_plus_rnaalifold_params)
  begin = time.time()
  pool.map(run_rnaalifold, mafft_xinsi_plus_rnaalifold_params)
  rnaalifold_elapsed_time = time.time() - begin
  pool.map(run_rnaalifold, ref_sa_plus_rnaalifold_params)
  print("The elapsed time of ConsAlifold = %f [s]." % consalifold_elapsed_time)
  print("The elapsed time of CentroidAlifold = %f [s]." % centroidalifold_elapsed_time)
  print("The elapsed time of RNAalifold = %f [s]." % rnaalifold_elapsed_time)
  print("The elapsed time of PETfold = %f [s]." % petfold_elapsed_time)

def run_mafft(mafft_params):
  (rna_seq_file_path, mafft_output_file_path) = mafft_params
  mafft_command = "mafft --thread 1 --quiet " + "--clustalout " + rna_seq_file_path + " > " + mafft_output_file_path
  utils.run_command(mafft_command)
  sa = AlignIO.read(mafft_output_file_path, "clustal")
  mafft_output_file_path = os.path.splitext(mafft_output_file_path)[0] + ".fa"
  AlignIO.write(sa, mafft_output_file_path, "fasta")

def run_probcons(probcons_params):
  (rna_seq_file_path, probcons_output_file_path) = probcons_params
  probcons_command = "probcons-RNA -clustalw " + rna_seq_file_path + " > " + probcons_output_file_path
  utils.run_command(probcons_command)
  sa = AlignIO.read(probcons_output_file_path, "clustal")
  AlignIO.write(sa, probcons_output_file_path, "clustal")
  probcons_output_file_path = os.path.splitext(probcons_output_file_path)[0] + ".fa"
  AlignIO.write(sa, probcons_output_file_path, "fasta")

def run_clustalw(clustalw_params):
  (rna_seq_file_path, clustalw_output_file_path) = clustalw_params
  clustalw_command = "clustalw -outorder=INPUT -infile=" + rna_seq_file_path + " -outfile=" + clustalw_output_file_path
  utils.run_command(clustalw_command)
  sa = AlignIO.read(clustalw_output_file_path, "clustal")
  clustalw_output_file_path = os.path.splitext(clustalw_output_file_path)[0] + ".fa"
  AlignIO.write(sa, clustalw_output_file_path, "fasta")

def run_mafft_xinsi(mafft_xinsi_params):
  (rna_seq_file_path, mafft_xinsi_output_file_path) = mafft_xinsi_params
  mafft_xinsi_command = "mafft-xinsi --thread 1 --quiet " + "--clustalout " + rna_seq_file_path + " > " + mafft_xinsi_output_file_path
  utils.run_command(mafft_xinsi_command)
  sa = AlignIO.read(mafft_xinsi_output_file_path, "clustal")
  mafft_xinsi_output_file_path = os.path.splitext(mafft_xinsi_output_file_path)[0] + ".fa"
  AlignIO.write(sa, mafft_xinsi_output_file_path, "fasta")

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

def run_petfold(petfold_params):
  (sa_file_path, petfold_output_file_path, gamma_str) = petfold_params
  petfold_command = "PETfold -f " + sa_file_path + " -a " + gamma_str
  (output, _, _) = utils.run_command(petfold_command)
  css = str(output).strip().split("\\n")[2].split("\\t")[1].strip()
  sta = AlignIO.read(sa_file_path, "fasta")
  AlignIO.write(sta, petfold_output_file_path, "stockholm")
  sta = AlignIO.read(petfold_output_file_path, "stockholm")
  sta.column_annotations["secondary_structure"] = css
  AlignIO.write(sta, petfold_output_file_path, "stockholm")

if __name__ == "__main__":
  main()
