***DON'T FORGET TO ADD THE WHOLE DIRECTORY AND SUBDIRECTORIES TO PATH IN MATLAB

=====================================================================================
FOLDERS DESCRIPTION:
=====================================================================================
DETECT_Versions - contains different versions of DETECT model
Plots - contains plotted outputs form, models
Plotting - contains code for easy 1 to 1 line plots
Additional Scripts - contain unrelated to the model but useful scripts.
=====================================================================================

=====================================================================================
MODEL VERESIONS DESCRIPTION (names are folder names at DETECT_Versions folder):
=====================================================================================
Original_DETECT - Original NSS, basic original model without changes

Efficient_DETECT - Efficient NSS, 60 times more efficient model version
*All models except Original are based on this version

FAKE_SS_DETECT - FAKE SS, keep SWC, Gness, and Tsoil constant

SS_DETECT - SS, Simplifed version with stable state solution

R_DETECT_NSS - Efficient NSS Rscript version
=====================================================================================

=====================================================================================
HOW TO RUN MODELS:
=====================================================================================
A. You can use actual code
1) Go to 'DETECT_Versions' folder *(cd comand in console or use matlab GUI) 
2) Go to 'xxx_DETECT' folder, where xxx is version of DETECT you interested in
*Note that R code version should be intiated differently.
3) Open 'xxx_ScriptFile.m' and Run it.
4) To check time took to run the model in console type 't2-t1'.
*Output will be saved to 'DETECT_Versions/xxx_DETECT/xxx_Outputs/output_6hourly.mat'
5) To plot output: from current folder go to 'xxx_Plots' folder, open, and...
... run 'xxx_Efficient_Rsoil_timeseries.m' file.

B. You can use 'Main_Program.m' script
1) Make Sure your current folder in MatLab is set to 'DETECT_MatLab'...
... and Run the 'Main_Program.m' script
2) In GUI choose model to run
3) In Parameter Values, choose parameters which you would like to use OR ignore.
4) Press 'Run Model' button and wait it to finish
*in MatLab windows you should be see it running(it finishes at 732 time steps)
*R version will not display anything in MatLab window, but will run in background
*time to run will be displayed in top right corner
5) To plot output: choose model you wish and press 'Plot Results'
*PLEASE NOTE THIS WAY TO RUN INGORCE ANY CHANGES IN 'ScriptFile.m' FILE TO...
YOU STILL CAN EDIT HARD CODE IN 'MainProgram.m'
=====================================================================================

=====================================================================================
DISCRIPTION OF MODEL FILES AND FOLDERS (original model used as example):
=====================================================================================
Folders:
Inputs - contain model inputs from PHACE site
Ouputs - contain model outputs
Plots - contain plots and plotting scripts

Files:
DETECT_Manual.docx - contains detailed information on DETECT model with how to run...
... it and adjust for different sites.
DETECTmodel.m - contains the main loop (the core) of DETECT model.
SourceTerm.m - contains submodels for microbial and root respiration.
ScriptFIle.m - contains initial parameters and options as well as save/load locations
=====================================================================================


=====================================================================================
TO MAKE SWC, T and Gness constant (original model used as example):
=====================================================================================
Go to Source_Term.m
After line 139:
fGness2=1+(Gness_IN-mean(Gness_IN))./(max(Gness_IN)-min(Gness_IN)); %used for v9g
ADD two lines:
fGness1 = mean(Gness_IN);
fGness2 = mean(Gness_IN);

Procced depending on how you run the model (ScriptFile or MainProgram):

A. If running through Script_File.m:
In Script_File.m
After line 55: SWC = Data.SWC;                %Nt * Nz
ADD: SWC(:,:) = mean(SWC(:));
After line 56: SoilT = Data.SoilT;            %Nt * Nz
ADD: SoilT(:,:) = mean(SoilT(:));

B. If running through MainProgram.m (currently only 4 models supported):

Depending from model find line which states:
Line 301: 'case 'Model_1','
- Original NSS, basic original model without changes

Line 322: 'case 'Model_2','
- Efficient NSS, 60 times more efficient model version

Line 343: 'case 'Model_3','
- FAKE SS, keep SWC, Gness, and Tsoil constant

Line 368: 'case 'Model_4','
- SS, Simplifed version with stable state solution

As you decide which model to run, after coresponding line add two lines:
SWC(:,:) = mean(SWC(:));
SoilT(:,:) = mean(SoilT(:));
=====================================================================================