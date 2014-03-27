#!/bin/sh

# CONFIG
DWI_BASE=/speechlab/subjects/pool_dwi
PROJ_PREFIX=RHY_
BATCH_DATA_BASE=/speechlab/5/scai/RHY/DATA_batch

SID=$1

FNIRT_DEST_DIR=${BATCH_DATA_BASE}/${SID}/fnirt
if [ ! -d ${FNIRT_DEST_DIR} ];
then
    echo Creating directory ${FNIRT_DEST_DIR}
    mkdir ${FNIRT_DEST_DIR}
fi

FNIRT_SRC_DIR=${DWI_BASE}/${PROJ_PREFIX}${SID}/fnirt
if [ ! -d ${FNIRT_SRC_DIR} ];
then 
    echo ERROR: Cannot find source fnirt directory: ${FNIRT_SRC_DIR}
    exit
fi


cp -v ${FNIRT_SRC_DIR}/fnirt_to_*.mat ${FNIRT_SRC_DIR}/fnirt_to_*.nii.gz ${FNIRT_SRC_DIR}/fnirt_to_*.msf ${FNIRT_DEST_DIR}