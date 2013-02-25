#!/usr/bin/python

import os
import sys
import glob
import argparse

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Put all files that match the input wild card into a list file")
    ap.add_argument("inWC", help="Input file pattern wild card")
    ap.add_argument("outListFN", help="Output list file name")

    if len(sys.argv) <= 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    inWC = args.inWC
    outListFN = args.outListFN

    d0 = glob.glob(inWC)
    d0.sort()

    outf = open(outListFN, "wt")
    for fn in d0:
        outf.write("%s\n" % fn)

    outf.close()
