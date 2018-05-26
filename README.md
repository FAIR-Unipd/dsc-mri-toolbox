
# Dynamic Susceptibility Contrast MRI toolbox

## About

Dynamic Susceptibility Contrast (DSC) MRI toolbox is a MATLAB toolbox to analyze DSC-MRI data. The code was implemented by [@marcocastellaro](https://github.com/marcocastellaro)  and Denis Peruzzo. This web page hosts the developmental source code. 

## Features:

### Pseudo-automatic AIF selection

The method is based on dicotomic hierarchical clustering method, please cite [1] if you use the AIF extraction feature:

[1] **Peruzzo Denis**,  Bertoldo Alessandra, Zanderigo Francesca and Cobelli Claudio, “[Automatic selection of arterial input function on dynamic contrast-enhanced MR images][paper1]”, *Computer methods and programs in biomedicine, 104:e148-e157 (2011)*.

### Deconvolution

  - SVD
  - cSVD
  - oSVD
  - NSR
  - Stable Spline
  
Please cite [2] if you perform the deconvolution with the Stable Spline algorithm.
  
[2] **Peruzzo Denis**,  Castellaro Marco, Pillonetto Gianluigi and Bertoldo Alessandra, “[Stable spline deconvolution for dynamic susceptibility contrast MRI][paper2]”, *Magnetic resonance in medicine, 78(5):1801--1811 (2017)*.
    
### Leakage correction


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

## License
This software is open source. The bulk of the code is covered by the MIT license. This tool is to be intended as a research tool and no medical decision should be made using it.



[//]: # (reference links)

   [paper1]: <https://www.sciencedirect.com/science/article/pii/S0169260711000447>

   [paper2]: <https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.26582>
