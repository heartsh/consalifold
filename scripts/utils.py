import os
import string
import subprocess
from Bio import AlignIO
from statistics import mode
import numpy

bracket_pairs = [("(", ")"), ("<", ">"), ("{", "}"), ("[", "]"), ("A", "a"), ("B", "b"), ("C", "c"), ("D", "d"), ("E", "e"), ]

def get_dir_paths():
  current_work_dir_path = os.getcwd()
  (head, tail) = os.path.split(current_work_dir_path)
  asset_dir_path = head + "/assets"
  program_dir_path = "/usr/local" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms"
  conda_program_dir_path = "/usr/local/ancnd/envs/rsrch" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms/ancnd/envs/rsrch"
  return (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path)

def get_css(css_file_path):
  sta = AlignIO.read(css_file_path, "stockholm")
  css_string = sta.column_annotations["secondary_structure"]
  sta_len = len(sta[0])
  num_of_rnas = len(sta)
  pos_map_sets = []
  for i in range(num_of_rnas):
    pos_map_sets.append([])
    pos = -1
    for j in range(sta_len):
      char = sta[i][j]
      if char != "-":
        pos += 1
      pos_map_sets[i].append(pos)
  css = []
  for i in range(num_of_rnas):
    css.append({})
  stack = []
  for (left, right) in bracket_pairs:
    for i, char in enumerate(css_string):
      if char == left:
        stack.append(i)
      elif char == right:
        col_pos = stack.pop()
        for j in range(num_of_rnas):
          base_pair_1 = (sta[j][col_pos], sta[j][i])
          if base_pair_1[0] == "-" or base_pair_1[1] == "-":
            continue
          pos_pair_1 = (pos_map_sets[j][col_pos], pos_map_sets[j][i])
          css[j][pos_pair_1] = True
  return css

def run_command(command):
  subproc = subprocess.Popen(
    command,
    stdout = subprocess.PIPE,
    shell = True
  )
  (output, error) = subproc.communicate()
  returned_code = subproc.wait()
  return (output, error, returned_code)

def get_avg_bpp_mat(bpp_mat_file_path, seq_lens, pos_map_sets, sta):
  sta_len = len(sta[0])
  bpp_mats = get_bpp_mats(bpp_mat_file_path, seq_lens)
  avg_bpp_mat = numpy.zeros((sta_len, sta_len))
  for i in range(sta_len):
    for j in range(i + 1, sta_len):
      avg_bpp = 0.
      effect_num_of_rnas = 0
      for k, pos_maps in enumerate(pos_map_sets):
        if sta[k][i] == "-" or sta[k][j] == "-":
          continue
        effect_num_of_rnas += 1
        bpp_mat = bpp_mats[k]
        pos_pair = (pos_maps[i], pos_maps[j])
        avg_bpp += bpp_mat[pos_pair]
      if effect_num_of_rnas > 0:
        avg_bpp /= effect_num_of_rnas
      avg_bpp_mat[i, j] = avg_bpp
  return avg_bpp_mat

def get_avg_upp_mat(upp_mat_file_path, seq_lens, pos_map_sets, sta):
  sta_len = len(sta[0])
  upp_mats = get_upp_mats(upp_mat_file_path, seq_lens)
  avg_upp_mat = numpy.zeros(sta_len)
  for i in range(sta_len):
    avg_upp = 0.
    effect_num_of_rnas = 0
    for k, pos_maps in enumerate(pos_map_sets):
      if sta[k][i] == "-":
        continue
      effect_num_of_rnas += 1
      upp_mat = upp_mats[k]
      pos = pos_maps[i]
      avg_upp += upp_mat[pos]
    if effect_num_of_rnas > 0:
      avg_upp /= effect_num_of_rnas
    avg_upp_mat[i] = avg_upp
  return avg_upp_mat

def get_bpp_mats(bpp_mat_file_path, seq_lens):
  bpp_mats = {}
  bpp_mat_file = open(bpp_mat_file_path)
  lines = bpp_mat_file.readlines()
  lines = [line for line in lines if line[0].isdigit() or line[0].startswith(">")]
  num_of_lines = len(lines)
  for i in range(0, num_of_lines - 1, 2):
    rna_id = int(lines[i][1 :])
    seq_len = seq_lens[rna_id]
    bpp_mat = numpy.zeros((seq_len, seq_len))
    for string in lines[i + 1].strip().split(" "):
      substrings = string.split(",")
      (j, k, bpp) = (int(substrings[0]), int(substrings[1]), max(0, min(1, float(substrings[2]))))
      bpp_mat[j, k] = bpp
    bpp_mats[rna_id] = bpp_mat
  return bpp_mats

def get_upp_mats(upp_mat_file_path, seq_lens):
  upp_mats = {}
  upp_mat_file = open(upp_mat_file_path)
  lines = upp_mat_file.readlines()
  lines = [line for line in lines if line[0].isdigit() or line[0].startswith(">")]
  num_of_lines = len(lines)
  for i in range(0, num_of_lines - 1, 2):
    rna_id = int(lines[i][1 :])
    seq_len = seq_lens[rna_id]
    upp_mat = numpy.zeros(seq_len)
    for string in lines[i + 1].strip().split(" "):
      substrings = string.split(",")
      (j, upp) = (int(substrings[0]), max(0, min(1, float(substrings[1]))))
      upp_mat[j] = upp
    upp_mats[rna_id] = upp_mat
  return upp_mats

def get_capr_prof_seqs(capr_output_file_path):
  capr_prof_seqs = {}
  capr_output_file = open(capr_output_file_path)
  lines = capr_output_file.readlines()
  lines = [line for line in lines]
  num_of_lines = len(lines)
  for line in lines[1:]:
    strings = line.split()
    prof_type = strings[0]
    profs = []
    for string in strings[1:]:
      profs.append(float(string))
    capr_prof_seqs[prof_type] = profs
  return capr_prof_seqs
