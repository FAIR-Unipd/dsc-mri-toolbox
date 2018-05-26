
# Dynamic Susceptibility Contrast MRI toolbox

## Features:

### Pseudo-automatic AIF selection

### Deconvolution

  - SVD
  - cSVD
  - oSVD
  - NSR
  - Stable Spline
### Leakage correction

This tool is to be intended as a research tool and no medical decision should be made using it.

Please cite if you use the Arterial Input Function tool 

```
@article{peruzzo2011automatic,
  title={Automatic selection of arterial input function on dynamic contrast-enhanced MR images},
  author={Peruzzo, Denis and Bertoldo, Alessandra and Zanderigo, Francesca and Cobelli, Claudio},
  journal={Computer methods and programs in biomedicine},
  volume={104},
  number={3},
  pages={e148--e157},
  year={2011},
  publisher={Elsevier}
}
```

If you use the Stable Spline algorithm 

```
article{peruzzo2017stable,
  title={Stable spline deconvolution for dynamic susceptibility contrast Mri},
  author={Peruzzo, Denis and Castellaro, Marco and Pillonetto, Gianluigi and Bertoldo, Alessandra},
  journal={Magnetic resonance in medicine},
  volume={78},
  number={5},
  pages={1801--1811},
  year={2017},
  publisher={Wiley Online Library}
}
```

## Example

You need to load the dataset with matlab. First of all it is needed to convert the DSC acquisition from DICOM to NIfTI. We do suggest to use dicom2niix that can be downloaded [here](https://github.com/neurolabusc/dcm2niix). 

Matlab provides routines to load Nifti file however an alternative could be to use the NIfTI toolbox that can be downloaded [here](https://it.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image). 

Below the example code to load the data and perform a DSC perfusion quantification. By default the code will produce CBV,CBF and MTT maps. CBV will be also corrected for leackage if present. CBF by default will be computed by SVD,cSVD and oSVD and the correspondent MTT maps will be produced.

```
% load the dataset to be analyzed 
DSC_info   = niftiinfo(fullfile('demo-data','GRE_DSC.nii.gz'));
DSC_volume = niftiread(DSC_info);

% Set minimum acquistion parameters 
TE = 0.025; % 25ms
TR = 1.55;  % 1.55s

% Perform quantification
[cbv,cbf,mtt,cbv_lc]=DSC_mri_core(DSC_volume,TE,TR);
```

A good way to perform a debug of the process is to load the default options with ```DSC_mri_get_options``` and change the default display properties to 3.

```
custom_options = DSC_mri_get_options();
custom_options.display = 3;
[cbv,cbf,mtt,cbv_lc]=DSC_mri_core(DSC_volume,TE,TR,custom_options);
```
