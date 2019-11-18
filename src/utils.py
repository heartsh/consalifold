import os
import string
import subprocess
from Bio import AlignIO

def get_dir_paths():
  current_work_dir_path = os.getcwd()
  (head, tail) = os.path.split(current_work_dir_path)
  asset_dir_path = head + "/assets"
  program_dir_path = "/usr/local" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms"
  conda_program_dir_path = "/usr/local/ancnd/envs/rsrch" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms/ancnd/envs/rsrch"
  return (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path)

if False:
  def print_css(rna_id_pos_triple_seqs, seqs):
    num_of_seqs = len(seqs)
    css_strings = [list([" "] * len(seqs[i])) for i in range(0, num_of_seqs)]
    alphabet = list(string.printable)
    for (i, rna_id_pos_triples) in enumerate(rna_id_pos_triple_seqs):
      for (rna_id, pos_1, pos_2) in rna_id_pos_triples:
        char = alphabet[i]
        css_strings[rna_id][pos_1] = char
        css_strings[rna_id][pos_2] = char
    for i in range(0, num_of_seqs):
      print(str(seqs[i]))
      print("".join(css_strings[i]))

  def get_rna_id_pos_triple_seqs(css_file_path):
    rna_id_pos_triple_seqs = []
    css_file = open(css_file_path)
    lines = css_file.readlines()
    lines = [line for line in lines if line[0].isdigit()]
    for line in lines:
      rna_id_pos_triples = []
      for substr in line.split(" "):
        subsubstrs = substr.split(":")
        pos_pair = subsubstrs[1].split(",")
        rna_id_pos_triples.insert(0, (int(subsubstrs[0]), int(pos_pair[0]), int(pos_pair[1])))
      rna_id_pos_triple_seqs.insert(0, rna_id_pos_triples)
    return rna_id_pos_triple_seqs

def get_css_and_flat_css(css_file_path):
  sta = AlignIO.read(css_file_path, "stockholm")
  css_string = sta.column_annotations["secondary_structure"]
  num_of_rnas = len(sta)
  sta_len = len(sta[0])
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
  flat_css = []
  for i in range(num_of_rnas):
    if False:
      for j in range(i + 1, num_of_rnas):
        css[(i, j)] = {}
    flat_css.append({})
    css.append({})
  stack = []
  for i, char in enumerate(css_string):
    if char == "(" or char == "<" or char == "[" or char == "{":
      stack.append(i)
    elif char == ")" or char == ">" or char == "]" or char == "}":
      col_pos = stack.pop()
      for j in range(num_of_rnas):
        base_pair_1 = (sta[j][col_pos], sta[j][i])
        if base_pair_1[0] == "-" or base_pair_1[1] == "-":
          continue
        pos_pair_1 = (pos_map_sets[j][col_pos], pos_map_sets[j][i])
        flat_css[j][pos_pair_1[0]] = True
        flat_css[j][pos_pair_1[1]] = True
        css[j][pos_pair_1] = True
        if False:
          for k in range(j + 1, num_of_rnas):
            base_pair_2 = (sta[k][col_pos], sta[k][i])
            pos_pair_2 = (pos_map_sets[k][col_pos], pos_map_sets[k][i])
            if base_pair_2[0] == "-" or base_pair_2[1] == "-":
              continue
            css[(j, k)][(pos_pair_1[0], pos_pair_1[1], pos_pair_2[0], pos_pair_2[1])] = True
  return css, flat_css

def run_command(command):
  subproc = subprocess.Popen(
    command,
    stdout = subprocess.PIPE,
    shell = True
  )
  (output, error) = subproc.communicate()
  returned_code = subproc.wait()
  return (output, error, returned_code)
