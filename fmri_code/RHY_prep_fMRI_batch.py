#!/usr/bin/env python

import os
import sys
import glob
import argparse

from scai_utils import info_log, error_log, cmd_stdout, \
                       check_bin_path, check_file, check_dir, \
                       saydo

machineSettings = {"ba3": {"dataDir": "/users/cais/RHY/DATA", 
                           "batchDataDir": "/users/cais/RHY/DATA_batch", 
                           "fsDataDir": "/users/cais/RHY/FSDATA", 
                           "batchCfg": "/users/cais/RHY/ANALYSIS/RHY_fmri_surf.cfg", 
                           "boldVolWC": "{subjID}_bold_*s???a???.nii.gz", 
                           "modelName": "fmri"}, \
                   "_default_": {"dataDir": "/speechlab/5/scai/RHY/DATA", 
                                 "batchDataDir": "/speechlab/5/scai/RHY/DATA_batch", 
                                 "fsDataDir": "/speechlab/5/scai/RHY/FSDATA", 
                                 "batchCfg": "/speechlab/5/scai/RHY/ANALYSIS/RHY_fmri_surf_BU.cfg", 
                                 "boldVolWC": "{subjID}_bold_*s???a???.nii.gz", 
                                 "modelName": "fmri"}}

HEMIS = ["lh", "rh"]

if __name__ == "__main__":
    ap = argparse.ArgumentParser("Prepare for batch analysis of RHY fMRI data")
    ap.add_argument("subjID", help="Subject ID, with the prefix MRI_ included (e.g., MRI_AWS_M01)")
    ap.add_argument("--run-batch", dest="bRunBatch", action="store_true",
                    help="Run the batch commands automatically (default: false)")
    
    if len(sys.argv) == 1:
        ap.print_help()
        sys.exit(0)

    args = ap.parse_args()

    machInfo = cmd_stdout("uname -a")
    hostName = machInfo[0].split(" ")[1]

    info_log("hostName = %s" % hostName)
    
    #=== Check the validity of subjID ===#
    if not args.subjID.startswith("MRI_"):
        info_log("subjID should usually start with the prefix MRI_. This is not the case in the entered subject ID (%s)" % (args.subjID), 
                 bWarn=True)
    
    sID = args.subjID.replace("MRI_", "")

    if machineSettings.keys().count(hostName) == 0:
        info_log("Cannot find host name %s in machineSettings. _default_ settings will be used" % hostName)
        hostName = "_default_"
        
    #=== Set up some paths ===#
    dataDir = machineSettings[hostName]["dataDir"]
    batchDataDir = machineSettings[hostName]["batchDataDir"]

    fsDataDir = machineSettings[hostName]["fsDataDir"]

    sDataDir = os.path.join(dataDir, sID)
    sBatchDataDir = os.path.join(batchDataDir, sID)
    
    check_dir(dataDir)
    check_dir(sDataDir)
    check_dir(sBatchDataDir, bCreate=True)
    check_dir(fsDataDir)

    boldDir = os.path.join(sDataDir, "bold")
    check_dir(boldDir)

    #=== Look for the fmri_model.mat and fmri_contrasts.mat files ===#
    modelMat = os.path.join(sDataDir, "fmri_model.mat")
    check_file(modelMat)
    info_log("modelMat = %s" % modelMat)

    contrMat = os.path.join(sDataDir, "fmri_contrasts.mat")
    check_file(contrMat)
    info_log("contrMat = %s" % contrMat)

    #=== Copy over the .mat files ===#
    saydo("cp %s %s/" % (modelMat, sBatchDataDir))
    saydo("cp %s %s/" % (contrMat, sBatchDataDir))

    #=== Load model mat and determine how many runs were collected ===#
    from scipy.io import loadmat
    
    mdl = loadmat(modelMat)
    nRuns = len(mdl["sess"][0])

    info_log("model file %s contains %d runs." % (modelMat, nRuns))

    #=== Look for raw bold volumes ===#
    info_log("Raw BOLD image dir = %s" % boldDir)        

    boldWC = machineSettings[hostName]["boldVolWC"].replace("{subjID}", sID)
    boldWC = os.path.join(boldDir, boldWC)
    dfs = glob.glob(boldWC)
    dfs.sort()

    if len(dfs) == 0:
        raise Exception, "Cannot find any BOLD volumes that match the path wild card: %s" % boldWC

    if len(dfs) == nRuns:
        info_log("The number of image files matches nRuns. Good")
        boldVols = dfs
    elif len(dfs) == 2 * nRuns:
        info_log("The number of image files equals 2 x nRuns. Assuming that thefirst of each pair os the non-motion corected volume.")
        boldVols = dfs[0::2]
    else:
        raise Exception, "The number of .nii.gz serie in directory %s is neither equal to nRuns=%d or 2*nRuns=%d" % (boldDir, nRuns, 2 * nRuns)

    #=== Copy image files and convert them to .nii ===#
    niiDir = os.path.join(sBatchDataDir, "nii")
    check_dir(niiDir, bCreate=True)

    check_bin_path("mri_convert")

    niftiFNs = []
    for i0 in range(nRuns):
        niftiFNs.append(os.path.join(niiDir, "00-%d.nii" % (i0 + 1)))
        saydo("mri_convert %s %s" % (boldVols[i0], niftiFNs[-1]))
        check_file(niftiFNs[-1])

    #=== Determine the number of contrasts ===#
    contr = loadmat(contrMat)

    nContrasts = len(contr["RegNameTContrasts"][0])
    info_log("Number of contrasts = %d" % nContrasts)


    #=== Check for the FS subject data dir ===#
    sFSDir = os.path.join(fsDataDir, sID)
    check_dir(sFSDir)

    envSubjsDir = os.getenv("SUBJECTS_DIR")
    if envSubjsDir != fsDataDir:
        info_log("The environmental variable SUBJECTS_DIR (=%s) does not match that in machineSettings (=%s). Set it properly before running the batch analysis" % (envSubjsDir, fsDataDir),                
                 bWarn=True)
        info_log("# Suggested command: ", bWarn=True)
        info_log("\texport SUBJECTS_DIR=%s" % os.path.abspath(fsDataDir), 
                 bWarn=True)

    #=== Create func_runs.txt file ===#
    scriptsDir = os.path.join(sBatchDataDir, "scripts")
    check_dir(scriptsDir, scriptsDir)

    funcRunsFN = os.path.join(scriptsDir, "func_runs.txt")
    funcRunsF = open(funcRunsFN, "wt")
    for i0 in range(nRuns):
        funcRunsF.write("%d " % (i0 + 1))
    funcRunsF.write("\n")
    funcRunsF.close()
    
    check_file(funcRunsFN)
    info_log("Run numbers written to file: %s" % funcRunsFN)

    #=== Check batch script prerequisites ===#
    batchCfg = os.path.abspath(machineSettings[hostName]["batchCfg"])
    check_file(batchCfg)

    batchSteps = ["motion_correct", "smooth", "preoutlier_model", 
                  "outliers", "prepare_analysis", "run_firstlevel", 
                  "run_contrasts"]
    batchCmds = []

    info_log("# To run the batch analysis, first do: ")
    info_log("\tcd %s" % (sBatchDataDir))

    for (i0, step) in enumerate(batchSteps):
        tCmd = "run_subject.sh -c %s -s %s -S %s -m %s %s" \
               % (batchCfg, step, step, \
                  machineSettings[hostName]["modelName"], sID)
        batchCmds.append(tCmd)
        info_log("# Command for %s: " % step)
        info_log("\t%s" % batchCmds[-1])

    #=== Prepare bbregister command ===#
    meanBoldVol = os.path.join(niiDir, "mean00-1.nii")
    func2Struct_dat = os.path.join(niiDir, "func2struct.bbr.dat")
    func2Struct_mat = os.path.join(niiDir, "func2struct.bbr.mat")
    bbrCmd = "bbregister --s %s --mov %s --reg %s --fslmat %s --bold --init-fsl" \
             % (sID, meanBoldVol, func2Struct_dat, func2Struct_mat)
    info_log("# Command for bbregister: ")
    info_log("\t%s" % bbrCmd)

    #=== Prepare contrast viewing command ===#
    L1Dir = os.path.join(sBatchDataDir, \
                         "firstlevel_%s" \
                         % machineSettings[hostName]["modelName"])


    for i0 in range(nContrasts):
        spmTViewCmd_vol = "tkmedit %s T1.mgz -surfs -overlay %s -overlay-reg %s -fthresh 6 -fmid 9" \
                          % (sID, os.path.join(L1Dir, "spmT_%.4d.img" % (i0 + 1)), \
                             func2Struct_dat)

        info_log("# Commands for viewing spmT for contrast #%d in the volume: " % (i0 + 1))
        info_log("\t%s" % spmTViewCmd_vol)

        info_log("# Command for viewing spmT for contrast #%d on the surface: " \
                 % (i0 + 1))
        
        for hemi in HEMIS:
            spmTViewCmd_surf = "tksurfer %s %s inflated -gray -overlay %s -ovelay-reg %s -fthresh 6 -fmid 9" \
                               % (sID, hemi,\
                                  os.path.join(L1Dir, \
                                               "spmT_%.4d.img" % (i0 + 1)),\
                                  func2Struct_dat)
            info_log("\t%s" % spmTViewCmd_surf)

        info_log(" ")

    #=== (Optional): Automatically run the batch commands ===%
    if args.bRunBatch:
        for (i0, cmd) in enumerate(batchCmds):
            saydo(cmd)
