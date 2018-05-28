
# Dynamic Susceptibility Contrast MRI toolbox

## About

Dynamic Susceptibility Contrast (DSC) MRI toolbox is a MATLAB toolbox to analyze DSC-MRI data. The code was implemented by [@peruzzod](https://github.com/peruzzod) and [@marcocastellaro](https://github.com/marcocastellaro). This web page hosts the developmental source code. 

## Features

### Semi-automatic AIF selection

The method is based on dicotomic hierarchical clustering method, it only need to select the slice where it will look for the best AIF. Please cite [1] if you use the AIF extraction tool:

[1] **Peruzzo Denis**,  Bertoldo Alessandra, Zanderigo Francesca and Cobelli Claudio, “[Automatic selection of arterial input function on dynamic contrast-enhanced MR images][paper1]”, *Computer methods and programs in biomedicine, 104:e148-e157 (2011)*.

### Deconvolution

  - SVD
  - cSVD
  - oSVD
  - NSR

Please cite [2] if you perform the deconvolution with the Nonlinear Stochastic Regularization algorithm.
  
[2] **Zanderigo Francesca**, and Bertoldo Alessandra and Pillonetto Gianluigi and Cobelli Claudio, “[Nonlinear stochastic regularization to characterize tissue residue function in bolus-tracking MRI: assessment and comparison with SVD, block-circulant SVD, and Tikhonov][paper2]”, *IEEE Transactions on Biomedical Engineering, 56(5):1287--1297 (2009)*.  
  
  - Stable Spline
  
Please cite [3] if you perform the deconvolution with the Stable Spline algorithm.
  
[3] **Peruzzo Denis**,  Castellaro Marco, Pillonetto Gianluigi and Bertoldo Alessandra, “[Stable spline deconvolution for dynamic susceptibility contrast MRI][paper3]”, *Magnetic resonance in medicine, 78(5):1801--1811 (2017)*.
    
### Leakage correction

Leakage correction is implemented following the approach proposed by Boxerman et al. in this paper, please cite [4] if ou use this correction.

[4] **Boxerman JL**, Schmainda KM and Weisskoff RM “[Relative cerebral blood volume maps corrected for contrast agent extravasation significantly correlate with glioma tumor grade, whereas uncorrected maps do not][paper4]”, *American Journal of Neuroradiology, 27(4):859--867 (2006)*.

## Example

You need to load the dataset with matlab. First of all it is needed to convert the DSC acquisition from DICOM to NIfTI. We do suggest to use *dicom2niix* that can be downloaded [here](https://github.com/neurolabusc/dcm2niix). 

Matlab provides routines to load Nifti file however an alternative could be to use the *NIfTI toolbox* that can be downloaded [here](https://it.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image). 

Below the example code to load the data and perform a DSC perfusion quantification. By default the code will produce *CBV*, *CBF* and *MTT* maps. *CBV* will be also corrected for leackage if present. *CBF* by default will be computed by *SVD*, *cSVD* and *oSVD* and the correspondent *MTT* maps will be produced.

```
% load the dataset to be analyzed 
DSC_info   = niftiinfo(fullfile('demo-data','GRE_DSC.nii.gz'));
DSC_volume = niftiread(DSC_info);

% Set minimum acquistion parameters 
TE = 0.025; % 25ms
TR = 1.55;  % 1.55s

% Perform quantification
[cbv,cbf,mtt,cbv_lc,ttp,mask,aif,conc,s0]=DSC_mri_core(DSC_volume,TE,TR);
```

A good way to perform a debug of the process is to load the default options with ```DSC_mri_getOptions``` and change the default ```display``` properties to ```3```. 

Values for display verbose options 
 - ```0``` no diplay at all
 - ```1``` notify only with text messages
 - ```2``` notify with representative images
 - ```3``` debug - notify many images, to be used only for tests, please avoid id when analyzing many subjects

```
custom_options = DSC_mri_getOptions();
custom_options.display = 3;
[cbv,cbf,mtt,cbv_lc,ttp,mask,aif,conc,s0]=DSC_mri_core(DSC_volume,TE,TR,custom_options);
```

A demo file that perform analysis of a sample subject (included in the demo-data folder) can be find in the [DSC_main_demo.m](DSC_main_demo.m) file.

## GUI 

Once the quantification has been performed it is possible to use a very simple GUI with the command ```DSC_mri_show_results``` to display maps computed with the toolbox and also to evaluate how each method selected performed in term of residue function calculation and reconvolution with the AIF. It is possible to select the voxel to be inspected.

```
DSC_mri_show_results(cbv_lc,cbf,mtt,ttp,mask,aif,conc,s0);
```

## License
This software is open source. The bulk of the code is covered by the MIT license. This tool is to be intended as a research tool and no medical decision should be made using it.

[//]: # (reference links)

   [paper1]: <https://www.sciencedirect.com/science/article/pii/S0169260711000447>
   [paper2]: <https://ieeexplore.ieee.org/abstract/document/4770181/>
   [paper3]: <https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.26582>
   [paper4]: <http://www.ajnr.org/content/27/4/859.short>
   
