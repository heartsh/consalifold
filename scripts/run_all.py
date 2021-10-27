#! /usr/bin/env python

import compile_rna_fams
import grid_search
import run_css_estimation_programs
import run_css_estimation_programs_2
import get_stats_of_css_estimation_programs
import get_stats_of_css_estimation_programs_2
import get_stats_of_css_estimation_programs_3
import get_stats_of_css_estimation_programs_4
import conduct_running_time_comparison
import conduct_running_time_comparison_2
import get_prob_dists

def main():
  compile_rna_fams.main()
  grid_search.main()
  run_css_estimation_programs.main()
  run_css_estimation_programs_2.main()
  get_stats_of_css_estimation_programs.main()
  get_stats_of_css_estimation_programs_2.main()
  get_stats_of_css_estimation_programs_3.main()
  get_stats_of_css_estimation_programs_4.main()
  conduct_running_time_comparison.main()
  conduct_running_time_comparison_2.main()
  get_prob_dists.main()

if __name__ == "__main__":
  main()
