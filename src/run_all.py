#! /usr/bin/env python

import sample_rna_fams
import run_css_estimation_programs
import get_stats_of_css_estimation_programs

def main():
  # sample_rna_fams.main() # If you resample data, different figures to the figures shown in the paper might be reproduced.
  run_css_estimation_programs.main()
  get_stats_of_css_estimation_programs.main()

if __name__ == "__main__":
  main()
