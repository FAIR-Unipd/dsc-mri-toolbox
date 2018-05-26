function DSC_mri_save_nifti_from_hdr(img,old_nii,file_nii)


nifti = old_nii;
nifti.hdr.dime.datatype = 16;
nifti.img = img;
save_untouch_nii(nifti,file_nii);