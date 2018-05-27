% DSC_mri_toolbox demo

% load the dataset to be analyzed 
DSC_info   = niftiinfo(fullfile('demo-data','GRE_DSC.nii.gz'));
DSC_volume = niftiread(DSC_info);

% Set minimum acquistion parameters 
TE = 0.025; % 25ms
TR = 1.55;  % 1.55s

% Perform quantification
[cbv,cbf,mtt,cbv_lc,ttp,mask,aif,conc,s0]=DSC_mri_core(DSC_volume,TE,TR);

% View Results
DSC_mri_show_results(cbv_lc,cbf,mtt,ttp,mask,aif,conc,s0);
