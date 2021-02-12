#! /usr/bin/env python

import getopt, sys
from Bio import AlignIO
from RNA import fold_compound

def main():
  try:
    opts, args = getopt.getopt(sys.argv[1:], "", [])
  except getopt.GetoptError as err:
    print(err)
    usage()
    sys.exit(2)
  input_file_path = args[0]
  sa = AlignIO.read(input_file_path, "clustal")
  rows = [str(row.seq) for row in sa]
  fc = fold_compound(rows)
  fc.pf()
  bpp_mat = fc.bpp()
  sparse_bpp_mat = {}
  sa_len = len(sa[0])
  bpp_mat_str = ""
  for i in range(0, sa_len):
    for j in range(i + 1, sa_len):
      bpp = bpp_mat[i][j]
      if bpp > 0:
        bpp_mat_str += "%d,%d,%e " % (i - 1, j - 1, bpp)
  print(bpp_mat_str)

if __name__ == "__main__":
  main()
