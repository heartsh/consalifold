#! /usr/bin/env python

import get_comm_struct
import compile_rna_fams
import run_css_estimation_programs
import get_stats_of_css_estimation_programs

def main():
  get_comm_struct.main()
  # compile_rna_fams.main()
  run_css_estimation_programs.main()
  get_stats_of_css_estimation_programs.main()

if __name__ == "__main__":
  main()
