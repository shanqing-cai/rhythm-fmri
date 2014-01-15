#!/bin/bash

RHY_BASE=/speechlab/5/scai/RHY

rsync -avz cais@ba3.mit.edu:~/RHY/FSDATA ${RHY_BASE}/
rsync -avz ${RHY_BASE}/FSDATA cais@ba3.mit.edu:~/RHY/

rsync -avz cais@ba3.mit.edu:~/RHY/DATA_batch ${RHY_BASE}/
rsync -avz ${RHY_BASE}/DATA_batch cais@ba3.mit.edu:~/RHY/

rsync -avz cais@ba3.mit.edu:~/RHY/L2_batch ${RHY_BASE}/
rsync -avz ${RHY_BASE}/L2_batch cais@ba3.mit.edu:~/RHY/
