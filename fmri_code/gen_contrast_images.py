#!/usr/bin/env python

import sys
import os
import argparse
import tempfile

CONVERT_BIN = "convert"

from scai_utils import check_file, check_dir, check_bin_path, saydo

from fmri_settings import L2_BATCH_BASE, commonSurfID, VOL_TEMPLATE_FS_ID, \
                          HEMIS, CON_IMG_DIR, \
                          TKMEDIT_VOL_BRIGHTNESS, TKMEDIT_VOL_CONTRAST, \
                          TKMEDIT_CBM_SLICES, TKMEDIT_BGTHAL_SLICES, \
                          CBM_IMG_SIZE, CBM_IMG_OFFSET, \
                          BGTHAL_IMG_SIZE, BGTHAL_IMG_OFFSET
                          

import tempfile

def prep_tcl_script_surf(tclScript, imgFNs):
    tclF = open(tclScript, "wt")
    tclF.write("save_tiff %s\n" % imgFNs[0])
    tclF.write("rotate_brain_y 180\n");
    tclF.write("redraw\n")
    tclF.write("save_tiff %s\n" % imgFNs[1])
    tclF.write("exit\n")
    tclF.close()    

def prep_tcl_script_vol(tclScript, imgBase):
    tclF = open(tclScript, "wt")
    tclF.write("SetVolumeBrightnessContrast 0 %f %f\n" % \
               (TKMEDIT_VOL_BRIGHTNESS, TKMEDIT_VOL_CONTRAST))

    imgFNs = {"CBM": [None] * len(TKMEDIT_CBM_SLICES),
              "BGTHAL": [None] * len(TKMEDIT_BGTHAL_SLICES)}

    for (h0, zone) in enumerate(imgFNs.keys()):
        if zone == "CBM":
            slice_ns = TKMEDIT_CBM_SLICES
        elif zone == "BGTHAL":
            slice_ns = TKMEDIT_BGTHAL_SLICES
        else:
            raise Exception

        for (i0, slcn) in enumerate(slice_ns):
            imgFNs[zone][i0] = "%s_%s_%.3d.tiff" % (imgBase, zone, slcn)
            tclF.write("SetSlice %d\n" % slcn)
            tclF.write("RedrawScreen\n")
            tclF.write("SaveTIFF %s\n" % imgFNs[zone][i0])

    tclF.write("exit\n")

    return imgFNs

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Generate 2nd-level contrast images (RHY project)")
    ap.add_argument("grp", help="Group: {ALL, AWS, ANS}")
    ap.add_argument("contr", help="Contrast number (e.g., 0001)")
    ap.add_argument("type", help="Type: {sf10 | vol | both}")
    ap.add_argument("fThresh", type=float, 
                    help="Functional activation (sig) threshold (e.g., 3)")
    ap.add_argument("fMid", type=float, 
                    help="Functional activation (sig) mid value (e.g., 5)")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(1)

    args = ap.parse_args()

    if args.type == "both":
        types = ["sf10", "vol"]
    else:
        types = [args.type]

    check_bin_path(CONVERT_BIN)

    check_dir(CON_IMG_DIR)

    l2Base = os.path.abspath(L2_BATCH_BASE)
    check_dir(l2Base)
    
    stitchRes = {}
    for t_type in types:
        conDir = os.path.join(l2Base, 
                              "%s_%s_%s" % (args.grp, args.contr, t_type))
        check_dir(conDir)

        stitchRes[t_type] = {}

        if t_type.startswith("sf"):
            #=== OSGM and BGC ===#
            if args.grp == "ALL":
                grpContrasts = ["bgc", "osgm"]
            else:
                grpContrasts = ["osgm"]

            for (k0, grpCon) in enumerate(grpContrasts):
                stitchRes[t_type][grpCon] = "%s_%s_%s_%s.tiff" % \
                                    (args.grp, args.contr, 
                                     t_type, grpCon)
                stitchRes[t_type][grpCon] = \
                    os.path.join(CON_IMG_DIR, stitchRes[t_type][grpCon])
                imgFNs = {"lat": [], "med": []}
                for (i0, hemi) in enumerate(HEMIS):
                    imgFNs["lat"].append(os.path.join(CON_IMG_DIR,  
                                               "%s_%s_%s_%s.%s.lat.tiff" % \
                                               (args.grp, args.contr, t_type, 
                                                grpCon, hemi)))
                    imgFNs["med"].append(os.path.join(CON_IMG_DIR, 
                                               "%s_%s_%s_%s.%s.med.tiff" % \
                                               (args.grp, args.contr, t_type, 
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
                    t_imgFNs = []
                    t_imgFNs.append(imgFNs["lat"][-1])
                    t_imgFNs.append(imgFNs["med"][-1])
                    prep_tcl_script_surf(tclScript, t_imgFNs)

                    # Put in tcl script information 
                    tksurferCommand += "-tcl %s " % tclScript

                    # Execute tksurfer 
                    saydo(tksurferCommand)
                    saydo("rm -rf %s" % tclScript)
                    for imgFN in t_imgFNs:
                        check_file(imgFN)

                # Stitch images together
                # Horizontal stitches
                horizStitch = [""] * len(imgFNs.keys())
                for (j0, side) in enumerate(imgFNs.keys()):
                    horizStitch[j0] = tempfile.mktemp() + ".tiff"
                    stitchCmd = "%s %s %s +append %s" % (CONVERT_BIN, 
                                                         imgFNs[side][0], 
                                                         imgFNs[side][1], 
                                                         horizStitch[j0])
                    saydo(stitchCmd)
                    check_file(horizStitch[j0])

                stitchCmd = "%s " % CONVERT_BIN
                for (j0, stfn) in enumerate(horizStitch):
                    stitchCmd += "%s " % stfn
                stitchCmd += "-append %s" % stitchRes[t_type][grpCon]
                saydo(stitchCmd)
                check_file(stitchRes[t_type][grpCon])

                # Clean up
                for (j0, side) in enumerate(horizStitch):
                    saydo("rm -f %s" % horizStitch[j0])
                for (j0, side) in enumerate(imgFNs.keys()):
                    for (j1, fn) in enumerate(imgFNs[side]):
                        saydo("rm -f %s" % imgFNs[side][j1])

        elif t_type == "vol":
            #=== OSGM and BGC ===#
            if args.grp == "ALL":
                grpContrasts = ["bgc", "osgm"]
            else:
                grpContrasts = ["osgm"]

            for (k0, grpCon) in enumerate(grpContrasts):
                stitchRes[t_type][grpCon] = "%s_%s_%s_%s.tiff" % \
                                    (args.grp, args.contr, 
                                     t_type, grpCon)
                stitchRes[t_type][grpCon] = \
                    os.path.join(CON_IMG_DIR, stitchRes[t_type][grpCon])

                if grpCon == "osgm":
                    sigImg = os.path.join(conDir, grpCon, 
                                          "osgm", "sig.mgh")
                elif grpCon == "bgc":
                    sigImg = os.path.join(conDir, grpCon, 
                                          "con_01", "sig.mgh")
                check_file(sigImg)

                tkmeditCommand = "tkmedit %s brain.mgz " % VOL_TEMPLATE_FS_ID + \
                                 "-overlay %s -fthresh %f -fmid %f " % \
                                 (sigImg, args.fThresh, args.fMid)

                # Prepare the tcl script
                tclScript = tempfile.mktemp() + ".tcl"

                imgBase = tempfile.mktemp()
                imgFNs = prep_tcl_script_vol(tclScript, imgBase)

                # Put in tcl script information 
                tkmeditCommand += "-tcl %s " % tclScript

                # Execute tkmedit
                saydo(tkmeditCommand)

                for zone in imgFNs.keys():
                    for imgFN in imgFNs[zone]:
                        check_file(imgFN)

                # Clean up tcl script
                imgFNs_crop = {}
                zoneStitches = {}            
                mergeZoneCmd = "%s " % CONVERT_BIN
                for (i0, zone) in enumerate(imgFNs.keys()):
                    imgFNs_crop[zone] = [None] * len(imgFNs[zone])
                    zoneStitches[zone] = "%s_%s_%s_%s.%s.tiff" % \
                                        (args.grp, args.contr, 
                                         t_type, grpCon, zone)
                    mergeCmd = "%s " % CONVERT_BIN

                    if zone == "CBM":
                        img_size = CBM_IMG_SIZE
                        img_offset = CBM_IMG_OFFSET
                    elif zone == "BGTHAL":
                        img_size = BGTHAL_IMG_SIZE
                        img_offset = BGTHAL_IMG_OFFSET
                    else:
                        raise Exception

                    for (i1, fn) in enumerate(imgFNs[zone]):
                        imgFNs_crop[zone][i1] = fn.replace(".tiff", "_crop.tiff")

                        cropCmd = "convert %s -crop %dx%d+%d+%d %s" % \
                                  (fn, img_size[0], img_size[1], \
                                       img_offset[0], img_offset[1], \
                                   imgFNs_crop[zone][i1])
                        saydo(cropCmd)
                        check_file(imgFNs_crop[zone][i1])

                        mergeCmd += "%s " % imgFNs_crop[zone][i1]

                    mergeCmd += " -append %s" % zoneStitches[zone]
                    saydo(mergeCmd)
                    check_file(zoneStitches[zone])

                    mergeZoneCmd += "%s " % zoneStitches[zone]

                    # Clean up
                    """

                    """

                mergeZoneCmd += " +append %s" % stitchRes[t_type][grpCon]
                saydo(mergeZoneCmd)
                check_file(stitchRes[t_type][grpCon])

                # Clean up
                for zone in imgFNs.keys():
                    for imgFN in imgFNs[zone]:
                        saydo("rm -f %s" % imgFN)

                for zone in imgFNs_crop.keys():
                    for imgFN in imgFNs_crop[zone]:
                        saydo("rm -f %s" % imgFN)

                for zone in zoneStitches.keys():
                    saydo("rm -f %s" % zoneStitches[zone])


        # Print the location of the resultant images
        print("=== Resultant images: ===")
        for (i0, grpCon) in enumerate(grpContrasts):
            print("%s:\t%s" % (grpCon, stitchRes[t_type][grpCon]))

    
    #=== (Optional: args.type == "both"): combine surf and vol images ===#
    if args.type == "both":
        surfVolStitch = {}

        for grpCon in stitchRes["vol"].keys():
            surfVolStitch[grpCon] = \
                os.path.join(CON_IMG_DIR,  
                             "%s_%s_surfVol_%s.tiff" % \
                             (args.grp, args.contr, grpCon))
            combineCmd = "%s %s %s +append %s" % \
                         (CONVERT_BIN, stitchRes["sf10"][grpCon], 
                          stitchRes["vol"][grpCon], surfVolStitch[grpCon])
            saydo(combineCmd)
            check_file(surfVolStitch[grpCon])
            
        print("=== Resultant combined surface-volume images ===")
        for (i0, grpCon) in enumerate(surfVolStitch.keys()):
            print("%s:\t%s" % (grpCon, surfVolStitch[grpCon]))
    
