# rhythm-fmri
Code for stimulus presentation and online ASR for the RHYTHM-FMRI project (fMRI study on the neural bases of the rhythm effect in stuttering)

Author: Shanqing Cai (shanqing.cai@gmail.com)

USAGE:

The main program for experiment running is runExperiment.m

Examples:

runExperiment(config_file_name);

runExperiment(config_file_name, 'twoScreens'); # Presents the stimulus window on a separate monitor to the left of the main one.


NOTE: This code by itself cannot run properly. It requires 

1) the Mex program Audapter to process audio input and output. Write to the author to request Audapter if needed. 

2) Some scripts in S.Cai's commonmcode repository. 

3) The Julius / Julian ASR engine (Windows version), in conjunction with the American English acoustic models supplied by VoxForge: julius-3.5.2-quickstart-windows. Freely available from 
http://www.voxforge.org/home/download

4) Some ASR-related code in the asr and asrcode directories require the CMU American English pronunciation dictionary:
www.speech.cs.cmu.edu/cgi-bin/cmudict


Acknowledgement:
This project is supported by NIH grant R01-DC0007683 (PI: Frank Guenther).


