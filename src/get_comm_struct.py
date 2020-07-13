#! /usr/bin/env python

import utils
from Bio import SeqIO
import numpy
import seaborn
from matplotlib import pyplot
import matplotlib
import os
import multiprocessing
import time
import datetime
import shutil
import community
import networkx
node_size = 150
edge_width = 2
label_font_size = 10
edge_label_font_size = 4 
pyplot.figure(figsize = (12, 12))
cmap = pyplot.cm.viridis

def main():
  (current_work_dir_path, asset_dir_path, program_dir_path, conda_program_dir_path) = utils.get_dir_paths()
  image_dir_path = asset_dir_path + "/images"
  if not os.path.exists(image_dir_path):
    os.mkdir(image_dir_path)
  seq_file_path = asset_dir_path + "/homologs_of_pri_miR_16_2.fa"
  seqs = [rec for rec in SeqIO.parse(seq_file_path, "fasta")]
  seq_lens = [len(seq) for seq in seqs]
  stat_dir_path = asset_dir_path + "/homologs_of_pri_miR_16_2"
  sta_file_path = stat_dir_path + "/gamma=8.sth"
  (_, _, css, flat_css, pos_map_sets, motifs, sta) = utils.get_css_and_flat_css(sta_file_path)
  sta_len = len(sta[0])
  access_bpp_mat_4_2l = utils.get_avg_bpp_mat(stat_dir_path + "/access_bpp_mats_on_2l.dat", seq_lens, pos_map_sets, sta)
  access_bpp_mat_4_ml = utils.get_avg_bpp_mat(stat_dir_path + "/access_bpp_mats_on_ml.dat", seq_lens, pos_map_sets, sta)
  bpp_mat_4_el = utils.get_avg_bpp_mat(stat_dir_path + "/bpp_mats_on_el.dat", seq_lens, pos_map_sets, sta)
  max_bpp_mat = numpy.zeros((sta_len, sta_len))
  for (i, j) in css.keys():
    max_bpp = 0
    bpp = access_bpp_mat_4_2l[i, j]
    if bpp > max_bpp:
      max_bpp = bpp
    bpp = access_bpp_mat_4_ml[i, j]
    if bpp > max_bpp:
      max_bpp = bpp
    bpp = bpp_mat_4_el[i, j]
    if bpp > max_bpp:
      max_bpp = bpp
    max_bpp_mat[i, j] = max_bpp
  graph = networkx.Graph()
  graph.add_nodes_from([i for i in range(sta_len + 5)])
  pos = networkx.circular_layout(graph)
  labels = {}
  edges = []
  edge_weights = []
  for (i, j) in css.keys():
    max_bpp = max_bpp_mat[i, j]
    if max_bpp == access_bpp_mat_4_2l[i, j]:
      edges.append((i, j))
      edge_weights.append(max_bpp)
  networkx.draw_networkx_edges(graph, pos, edgelist = edges, edge_color = edge_weights, edge_cmap = cmap, style = "solid", width = edge_width, edge_vmin = 0, edge_vmax = 1)
  edges = []
  edge_weights = []
  for (i, j) in css:
    max_bpp = max_bpp_mat[i, j]
    if max_bpp == access_bpp_mat_4_ml[i, j]:
      edges.append((i, j))
      edge_weights.append(max_bpp)
  networkx.draw_networkx_edges(graph, pos, edgelist = edges, edge_color = edge_weights, edge_cmap = cmap, style = "dashed", width = edge_width, edge_vmin = 0, edge_vmax = 1)
  edges = []
  edge_weights = []
  for (i, j) in css:
    max_bpp = max_bpp_mat[i, j]
    if max_bpp == bpp_mat_4_el[i, j]:
      edges.append((i, j))
      edge_weights.append(max_bpp)
  networkx.draw_networkx_edges(graph, pos, edgelist = edges, edge_color = edge_weights, edge_cmap = cmap, style = "dotted", width = edge_width, edge_vmin = 0, edge_vmax = 1)
  pyplot.colorbar(pyplot.cm.ScalarMappable(norm = matplotlib.colors.Normalize(), cmap = cmap), fraction = 0.01)
  upp_mat_4_2l = utils.get_avg_upp_mat(stat_dir_path + "/upp_mats_on_2l.dat", seq_lens, pos_map_sets, sta)
  upp_mat_4_ml = utils.get_avg_upp_mat(stat_dir_path + "/upp_mats_on_ml.dat", seq_lens, pos_map_sets, sta)
  upp_mat_4_el = utils.get_avg_upp_mat(stat_dir_path + "/upp_mats_on_el.dat", seq_lens, pos_map_sets, sta)
  upp_mat_4_hl = utils.get_avg_upp_mat(stat_dir_path + "/upp_mats_on_hl.dat", seq_lens, pos_map_sets, sta)
  max_upp_mat = numpy.zeros(sta_len)
  for i in range(sta_len):
    max_upp = 0
    upp = upp_mat_4_hl[i]
    if upp > max_upp:
      max_upp = upp
    upp = upp_mat_4_2l[i]
    if upp > max_upp:
      max_upp = upp
    upp = upp_mat_4_ml[i]
    if upp > max_upp:
      max_upp = upp
    upp = upp_mat_4_el[i]
    if upp > max_upp:
      max_upp = upp
    max_upp_mat[i] = max_upp
  nodes = []
  node_weights = []
  for i in range(sta_len):
    labels[i] = motifs[i]
    max_upp = max_upp_mat[i]
    if (not i in flat_css) and max_upp == upp_mat_4_2l[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "o", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  nodes = []
  node_weights = []
  for i in range(sta_len):
    max_upp = max_upp_mat[i]
    if (not i in flat_css) and max_upp == upp_mat_4_ml[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "p", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  nodes = []
  node_weights = []
  for i in range(sta_len):
    max_upp = max_upp_mat[i]
    if (not i in flat_css) and max_upp == upp_mat_4_el[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "v", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  nodes = []
  node_weights = []
  for i in range(sta_len):
    max_upp = max_upp_mat[i]
    if(not i in flat_css) and  max_upp == upp_mat_4_hl[i]:
      nodes.append(i)
      node_weights.append(max_upp)
  networkx.draw_networkx_nodes(graph, pos, nodelist = nodes, node_color = node_weights, node_shape = "d", cmap = cmap, node_size = node_size, vmin = 0, vmax = 1)
  networkx.draw_networkx_labels(graph, pos, labels = labels, font_color = "r", font_size = label_font_size)
  pyplot.axis("off")
  pyplot.savefig(image_dir_path + "/comm_struct.eps", bbox_inches = "tight")

if __name__ == "__main__":
  main()
