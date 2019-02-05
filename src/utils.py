import os
import string

def get_dir_paths():
  current_work_dir_path = os.getcwd()
  (head, tail) = os.path.split(current_work_dir_path)
  asset_dir_path = head + "/assets"
  program_dir_path = "/usr/local" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms"
  conda_program_dir_path = "/usr/local/ancnd/envs/rsrch" if current_work_dir_path.find("/home/masaki") == -1 else "/home/masaki/prgrms/ancnd/envs/rsrch"
  return (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path)

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
