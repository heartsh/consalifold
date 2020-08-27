#! /usr/bin/env python

import utils
from Bio import SeqIO
import numpy

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  input_rna_file_path = asset_dir_path + "/trna.fa"
  num_of_samples = 6
  records = numpy.random.choice([record for record in SeqIO.parse(input_rna_file_path, "fasta")], num_of_samples, replace = False)
  output_rna_file_path = asset_dir_path + "/sampled_trnas.fa"
  SeqIO.write(records, output_rna_file_path, "fasta")

if __name__ == "__main__":
  main()
