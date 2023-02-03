#! /usr/bin/env sh

./compile_rna_fams.py \
  && ./ grid_search.py \
  && ./run_css_estimation_programs.py \
  && ./run_css_estimation_programs_2.py \
  && ./get_stats_of_css_estimation_programs.py \
  && ./get_stats_of_css_estimation_programs_2.py \
  && ./get_stats_of_css_estimation_programs_3.py \
  && ./conduct_running_time_comparison.py \
  && ./conduct_running_time_comparison_2.py \
  && ./get_prob_dists.py
