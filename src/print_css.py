#! /usr/bin/env python

import utils
from Bio import SeqIO

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  rna_file_path = asset_dir_path + "/sampled_trnas.fa"
  css_file_path = asset_dir_path + "/sampled_trnas.dat"
  seqs = [record.seq for record in SeqIO.parse(rna_file_path, "fasta")]
  rna_id_pos_triple_seqs = utils.get_rna_id_pos_triple_seqs(css_file_path)
  utils.print_css(rna_id_pos_triple_seqs, seqs)

if __name__ == "__main__":
  main()
