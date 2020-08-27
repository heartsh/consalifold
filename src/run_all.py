#! /usr/bin/env python

import compile_rna_fams
import run_css_estimation_programs
import get_stats_of_css_estimation_programs
import sample_rnas
import get_trna_prob_dists

def main():
  compile_rna_fams.main()
  run_css_estimation_programs.main()
  get_stats_of_css_estimation_programs.main()
  # sample_rnas.main()
  get_trna_prob_dists.main()

if __name__ == "__main__":
  main()
