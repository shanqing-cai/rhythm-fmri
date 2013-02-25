#!/usr/bin/python

import os
import sys
import argparse

BASE_JCONF = "/home/cais/asr/scripts/julian.jconf"
HMMDEF = "/home/cais/asr/julius-3.5.2-quickstart-linux/acoustic_model_files_build726/hmmdefs"
TIEDLIST = "/home/cais/asr/julius-3.5.2-quickstart-linux/acoustic_model_files_build726/tiedlist"

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Generate the jconf file based on a based jconf file")
    ap.add_argument("inDir", help="Input directory")
    
    if len(sys.argv) <= 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()
    inDir = args.inDir

    # Copy the tiedlist to target directory
    os.system("cp %s %s" % (TIEDLIST, inDir))
    (tlpath, tlfn) = os.path.split(TIEDLIST)
    tiedList = os.path.join(inDir, tlfn)

    # Copy the base jconf to the target directory
    os.system("cp %s %s" % (BASE_JCONF, inDir))
    (jcpath, jcfn) = os.path.split(BASE_JCONF)
    jconf = os.path.join(inDir, jcfn)

    # Read and modify jconf text
    jconff = open(jconf, "rt")
    jct = jconff.read().split('\n')
    jconff.close()

    for (i0, t_line) in enumerate(jct):
        if t_line.startswith("-h "):
            jct[i0] = "-h %s" % HMMDEF
        elif t_line.startswith("-hlist "):
            jct[i0] = "-hlist %s" % tiedList
        elif t_line.startswith("-dfa "):
            jct[i0] = "-dfa %s" % os.path.join(inDir, "gram.dfa")
        elif t_line.startswith("-v "):
            jct[i0] = "-v %s" % os.path.join(inDir, "gram.dict")
            
    out_jconf = os.path.join(inDir, "jconf")
    out_jconf_f = open(out_jconf, "wt")
    for (i0, t_line) in enumerate(jct):
        out_jconf_f.write("%s\n" % t_line)
    out_jconf_f.close()
    
