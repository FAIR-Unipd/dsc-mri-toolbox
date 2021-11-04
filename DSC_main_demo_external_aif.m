% DSC_mri_toolbox demo

% ------ Load the dataset to be analyzed ---------------------------------
DSC_info   = niftiinfo(fullfile('demo-data','GRE_DSC.nii.gz'));
DSC_volume = niftiread(DSC_info);

% ------ Load pre-calculated AIF -----------------------------------------
load('demo-data/aif_conc.mat','aif_conc')

% ------ Set minimum acquistion parameters -------------------------------
TE = 0.025; % 25ms
TR = 1.55;  % 1.55s

% ------ Perform quantification ------------------------------------------ 
% Input   DSC_volume (4D matrix with raw GRE-DSC acquisition)
%         TE         (Echo time)
%         TR         (Repetition time)
% Output  cbv        (3D matrix with standard rCBV values)
%         cbf        (struct with 3D matrices of rCBF values for each method selected)
%         mtt        (struct with 3D matrices of MTT values for each method selected)
%         cbv_lc     (3D matrix with leackage corrected rCBV values)
%         ttp        (3D matrix with leackage corrected Time to Peak values)
%         mask       (3D matrix with computed mask)
%         aif        (struct with AIF extracted with clustering algorithm)
%         conc       (4D matrix with pseudo-concentration values)
%         s0         (3D matrix with S0 estimates from pre-contrast images)

opt_dsc_MRI=DSC_mri_getOptions();
opt_dsc_MRI.aif.enable= 0; 
aif.fit.gv = aif_conc;

[cbv,cbf,mtt,cbv_lc,ttp,mask,aif,conc,s0]=DSC_mri_core(DSC_volume,TE,TR,opt_dsc_MRI,aif);

% ------  View Results --------------------------------------------------- 
DSC_mri_show_results(cbv_lc,cbf,mtt,ttp,mask,aif,conc,s0);
