#!/usr/bin/env python

import sys
import os
import argparse
import tempfile

from scai_utils import check_file, check_dir, saydo

from fmri_settings import L2_BATCH_BASE, commonSurfID, VOL_TEMPLATE_FS_ID, \
                          HEMIS, CON_IMG_DIR

def prep_tcl_script(tclScript, imgFNs):
    tclF = open(tclScript, "wt")
    tclF.write("save_tiff %s\n" % imgFNs[0])
    tclF.write("rotate_brain_y 180\n");
    tclF.write("redraw\n")
    tclF.write("save_tiff %s\n" % imgFNs[1])
    tclF.write("exit\n")
    tclF.close()    

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Generate 2nd-level contrast images (RHY project)")
    ap.add_argument("grp", help="Group: {ALL, AWS, ANS}")
    ap.add_argument("contr", help="Contrast number (e.g., 0001)")
    ap.add_argument("type", help="Type: {sf10 | vol}")
    ap.add_argument("fThresh", type=float, 
                    help="Functional activation (sig) threshold (e.g., 3)")
    ap.add_argument("fMid", type=float, 
                    help="Functional activation (sig) mid value (e.g., 5)")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(1)

    args = ap.parse_args()

    check_dir(CON_IMG_DIR)

    l2Base = os.path.abspath(L2_BATCH_BASE)
    check_dir(l2Base)
    
    conDir = os.path.join(l2Base, 
                          "%s_%s_%s" % (args.grp, args.contr, args.type))
    check_dir(conDir)

    if args.type.startswith("sf"):
        #=== OSGM and BGC ===#
        if args.grp == "ALL":
            grpContrasts = ["bgc", "osgm"]
        else:
            grpContrasts = ["osgm"]

        for (k0, grpCon) in enumerate(grpContrasts):
            for (i0, hemi) in enumerate(HEMIS):
                imgFNs = []
                imgFNs.append(os.path.join(CON_IMG_DIR,  
                                           "%s_%s_%s_%s.%s.lat.tiff" % \
                                           (args.grp, args.contr, args.type, 
                                            grpCon, hemi)))
                imgFNs.append(os.path.join(CON_IMG_DIR, 
                                           "%s_%s_%s_%s.%s.med.tiff" % \
                                           (args.grp, args.contr, args.type, 
                                            grpCon, hemi)))

                if grpCon == "osgm":
                    sigImg = os.path.join(conDir, "%s_%s" % (grpCon, hemi), 
                                          grpCon, "sig.mgh")
                elif grpCon == "bgc":
                    sigImg = os.path.join(conDir, "%s_%s" % (grpCon, hemi), 
                                          "con_01", "sig.mgh")
                check_file(sigImg)

                tksurferCommand = "tksurfer %s %s inflated -gray " % \
                                  (commonSurfID, hemi, ) + \
                                  "-overlay %s -fthresh %f -fmid %f " % \
                                  (sigImg, args.fThresh, args.fMid) + \
                                  "-colscalebarflag 1 "

                # Prepare the tcl script
                tclScript = tempfile.mktemp() + ".tcl"
                prep_tcl_script(tclScript, imgFNs)
            
                # Put in tcl script information 
                tksurferCommand += "-tcl %s " % tclScript

                # Execute tksurfer 
                saydo(tksurferCommand)
                saydo("rm -rf %s" % tclScript)
                for imgFN in imgFNs:
                    check_file(imgFN)

    elif args.type == "vol":
        pass
    
