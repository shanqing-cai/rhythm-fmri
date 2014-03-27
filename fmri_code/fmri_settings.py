modelName = "fmri"
HEMIS = ["lh", "rh"]
commonSurfID = "fsaverage"

MATLAB_BIN = "matlab"

con_fns = {}
con_fns['10'] = 'con_10'
con_fns['01'] = 'con_01'
con_fns['0m1'] = '/users/cais/STUT/analysis/design/con_0-1'

VOL_TEMPLATE="/usr/share/fsl/5.0/data/standard/MNI152_T1_1mm_brain.nii.gz"
VOL_TEMPLATE_FS_ID = "MNI152_1mm"

BATCH_BASE = "/users/cais/RHY/DATA_batch"
L2_BATCH_BASE = "/users/cais/RHY/L2_batch"

CON_IMG_DIR = "/users/cais/RHY/con_images"

### Visualization options ###
# tkmedit brightness and contrast
TKMEDIT_VOL_BRIGHTNESS = 0.6
TKMEDIT_VOL_CONTRAST = 12

# Slice numbers for cerebellum (CBM) and BG-thalamus (BGTHAL)
TKMEDIT_CBM_SLICES = [60, 65, 70, 75, 80, 85, 90, 100, 105, 110]
TKMEDIT_BGTHAL_SLICES = [115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165]

CBM_IMG_SIZE = [300, 140]
CBM_IMG_OFFSET = [106, 300]

BGTHAL_IMG_SIZE = [200, 100]
BGTHAL_IMG_OFFSET = [156, 240]


