#!/usr/bin/env python

import os
import sys
import argparse
import glob

from scai_utils import check_file, check_dir, info_log, error_log, saydo

modelName = "fmri"
HEMIS = ["lh", "rh"]
commonSurfID = "fsaverage"

MATLAB_BIN = "matlab"

con_fns = {}
con_fns['10'] = 'con_10'
con_fns['01'] = 'con_01'
con_fns['0m1'] = '/users/cais/STUT/analysis/design/con_0-1'

def generate_design_matrix(isAWS, matFN, bReverse=False):
# Input arguments:
#     isAWS: vector of 0's and 1's
#     matFN: output design matrix .mat file name 
#     bReverse: if AWS should be set to -1 and ANS to +1 (default: no)

    if len(isAWS) < 2:
        error_log("The design matrix cannot be generated because there are fewer than two subjects")

    if not matFN.endswith(".mat"):
        error_log("The input .mat file name %s has a wrong extension name" % matFN)

    X_line = 'X = ['
    for (i0, t_isAWS) in enumerate(isAWS):
        t_x = (float(t_isAWS) * 2.0 - 1.0)
        if bReverse:
            X_line += "%.1f; " % -t_x
        else:
            X_line += "%.1f; " % t_x

    X_line = X_line[:-2] + "];\n"
    X_line += "X = [ones(%d , 1), X];\n" % len(isAWS)
    X_line += "save('%s', 'X', '-v4');\n" % os.path.abspath(matFN)
    print(X_line) # DEBUG
    
    (t_path, t_fn) = os.path.split(os.path.abspath(matFN))
    mScriptGenX = os.path.join(t_path, 'gen_%s' % t_fn.replace(".mat", ".m"))

    mScriptGenX_f = open(mScriptGenX, "wt")
    mScriptGenX_f.write(X_line)
    mScriptGenX_f.close()
    check_file(mScriptGenX)


    matlabCmd = "%s -nosplash -nodesktop -r 'run %s; exit; '" % \
                (MATLAB_BIN, mScriptGenX)
    saydo(matlabCmd)
    check_file(matFN)


if __name__ == "__main__":
    ap = argparse.ArgumentParser("2nd-level (group-level) fMRI analysis")
    ap.add_argument("contrastNum", type=int, help="Contrast number (e.g., 1, 2)")
    ap.add_argument("batchBase", help="Batch fMRI analysis base directory (e.g., ~/RHY/DATA_batch)")
    ap.add_argument("fsDir", help="FreeSurfer base directory: SUBJECTS_DIR (e.g., ~/RHY/FSDATA")
    ap.add_argument("resBase", help="Base directory for saving L2 analysis results")
    ap.add_argument("--group", dest="group", type=str, 
                    help="Subject group (e.g., AWS). If not specified, will include all subjects and perform between-group analysis")
    ap.add_argument("--fwhm-surf", dest="fwhmSurf", type=float, default=10, 
                    help="FWHM for surface-based smoothing (default=0)")
    ap.add_argument("--redo", dest="bRedo", action="store_true", 
                    help="Force redo all steps (default: false)")

    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)
        
    args = ap.parse_args()

    #=== Input sanity check ===#
    check_dir(args.batchBase)
    check_dir(args.fsDir)

    #=== Check that SUBJECTS_DIR matches input argument ===#
    envFSDir = os.path.abspath(os.getenv("SUBJECTS_DIR"))
    if envFSDir != os.path.abspath(args.fsDir):
        error_log("Input FreeSurfer SUBJECTS_DIR does not match environmental SUBECTS_DIR")

    #=== Get the list of subjects ===#
    if args.group == None or len(args.group) == 0:
        wc = "*"
        bAll = True
    else:
        wc = args.group + "*"
        bAll = False;

    ds = glob.glob(os.path.join(args.batchBase, wc))
    if len(ds) == 0:
        raise Exception, "Cannot find subjects that match the selection criterion in directory: %s" % args.batchBase

    subjIDs = []
    subjIsAWS = []
    ds.sort()
    for (i0, dn) in enumerate(ds):
        subjIDs.append(os.path.split(dn)[1])
        
        #== Check the existence of the FreeSurfer directory ==#
        check_dir(os.path.join(args.fsDir, subjIDs[-1]))

        if bAll:
            if subjIDs[-1].startswith("AWS_"):
                subjIsAWS.append(1)
            elif subjIDs[-1].startswith("ANS_"):
                subjIsAWS.append(0)
            else:
                error_log("Unrecognized subject ID prefix in %s" % subjIDs[-1])

    info_log("Found %d subjects" % len(subjIDs))
    if bAll:
        info_log("Number of AWS = %d" % subjIsAWS.count(1))
        info_log("Number of ANS = %d" % subjIsAWS.count(0))

        # Determine if between-group comparison is to be performed
        bGrpComp = not (subjIsAWS.count(1) == 0 or subjIsAWS.count(0) == 0)

    #=== Project the con values to the surface ===#
    commSurfConFNs = {}
    for hemi in HEMIS:
        commSurfConFNs[hemi] = []

    for sID in subjIDs:
        #== Locate the con file ==#
        L1Dir = os.path.join(args.batchBase, sID, "firstlevel_%s" % modelName)
        check_dir(L1Dir)
        
        conFN = os.path.join(L1Dir, "con_%.4d.img" % args.contrastNum)
        check_file(conFN)

        #=== Locate the bbregister file ===#
        bbrDat = os.path.join(args.batchBase, sID, "nii", "func2struct.bbr.dat")
        check_file(bbrDat)

        #== Use vol2surf command to project the contrast values to the subject's own surface ===#
        #== Then use mris_preproc to project the data from the subject's own surface to the common space (fsaverage) ==#
        for hemi in HEMIS:
            surfConFN = os.path.join(L1Dir, "con_%.4d.%s.mgz" % \
                                            (args.contrastNum, hemi))
            if args.bRedo or not os.path.isfile(surfConFN):
                projCmd = "mri_vol2surf --mov %s --reg %s --o %s --hemi %s --trgsubject %s --noreshape --interp trilin" % \
                          (conFN, bbrDat, surfConFN, hemi, sID)
                
                saydo(projCmd)
                check_file(surfConFN)
                
            commSurfConFN = os.path.join(L1Dir, "con_%.4d.%s.fsav.mgz" % \
                                                (args.contrastNum, hemi))
            if args.bRedo or not os.path.isfile(commSurfConFN):
                preprocCmd = "mris_preproc --s %s --hemi %s --is %s --fwhm %f --target %s --out %s" % \
                             (sID, hemi, surfConFN, args.fwhmSurf, \
                              commonSurfID, commSurfConFN)
                saydo(preprocCmd)
                check_file(commSurfConFN)
            
            commSurfConFNs[hemi].append(commSurfConFN)

    #=== Concatenate the files ===#
    # "sf" stands for surface FWHM
    if args.group != None and len(args.group) > 0:
        t_grp = args.group
    else:
        t_grp = "ALL"
    
    # Prepare output base directory
    outBase = os.path.join(args.resBase, "%s_%.4d_sf%d" % \
                           (t_grp, args.contrastNum, args.fwhmSurf))
    
    check_dir(outBase, bCreate=True)
    
    mergedSurfCons = {}
    for hemi in HEMIS:
        mergedSurfCon = os.path.join(outBase, "merged_con.%s.fsav.mgz" % hemi) 
        
        concatCmd = "mris_preproc --hemi %s --out %s --target %s " \
                    % (hemi, mergedSurfCon, commonSurfID)
        if args.bRedo or not os.path.isfile(mergedSurfCon):
            for (i0, csc) in enumerate(commSurfConFNs[hemi]):
                concatCmd += "--s %s --is %s " % (commonSurfID, csc)

            saydo(concatCmd)
            check_file(mergedSurfCon)

        mergedSurfCons[hemi] = mergedSurfCon

    #=== Prepare the group-level design matrix and contrast vectors ===#
    if bAll and bGrpComp:
        designDir = os.path.join(outBase, "design")
        check_dir(designDir, bCreate=True)
        
        X_mat = os.path.join(designDir, "X.mat")
        X_reverse_mat = os.path.join(designDir, "X_reverse.mat")

        if args.bRedo or not os.path.isfile(X_mat):
            generate_design_matrix(subjIsAWS, X_mat)
        if args.bRedo or not os.path.isfile(X_reverse_mat):
            generate_design_matrix(subjIsAWS, X_reverse_mat, bReverse=True)
    

    #=== Perform mri_glmfit (OSGM) ===#
    viewCmds_osgm = {}
    for hemi in HEMIS:
        glmDir = os.path.join(outBase, "osgm_%s" % hemi)
        sig_mgh = os.path.join(glmDir, "osgm", "sig.mgh")
        
        if args.bRedo or not os.path.isfile(sig_mgh):
            fitCmd = 'mri_glmfit --cortex --y %s --osgm --glmdir %s'\
                     % (mergedSurfCons[hemi], glmDir)
            saydo(fitCmd)
            check_file(sig_mgh)
        
        viewCmds_osgm[hemi] = \
            "tksurfer %s %s inflated -gray -overlay %s -fthresh 2.5" % \
            (commonSurfID, hemi, sig_mgh)

    #=== Perform mri_glmfit: between-group comparison ===#
    X_mats = [X_mat, X_reverse_mat]
    desNames = ["bgc", "bgc_rev"]
    viewCmds_bgc = {}
    for hemi in HEMIS:
        viewCmds_bgc[hemi] = {}
        for (i0, des_mat) in enumerate(X_mats):
            glmDir = os.path.join(outBase, "%s_%s" % (desNames[i0], hemi))
            check_file(con_fns["01"])
            sig_mgh = os.path.join(glmDir, con_fns["01"], "sig.mgh")
            
            if args.bRedo or not os.path.isfile(sig_mgh):
                fitCmd = 'mri_glmfit --cortex --y %s --X %s --C %s --glmdir %s'\
                     % (mergedSurfCons[hemi], des_mat, con_fns["01"], glmDir)
                
                saydo(fitCmd)
                check_file(sig_mgh)
                
        
            viewCmds_bgc[hemi][desNames[i0]] = \
                "tksurfer %s %s inflated -gray -overlay %s -fthresh 2.5" % \
                (commonSurfID, hemi, sig_mgh)
    
        
    #=== Print view commands: OSGM ===#
    print("\n\n=== Commands for viewing OSGM results ===")
    for hemi in HEMIS:
        print("Hemisphere %s: " % hemi)
        print("\t%s\n" % viewCmds_osgm[hemi])

    #=== Print view commands: between-group comparisons ===#
    if bAll and bGrpComp:
        print("\n\n=== Commands for viewing between-group comparison results ===")
        for hemi in HEMIS:
            for (i0, t_desName) in enumerate(desNames):
                print("Hemisphere %s - %s: " % (hemi, t_desName))
                print("\t%s\n" % viewCmds_bgc[hemi][t_desName])
        print("\n")

