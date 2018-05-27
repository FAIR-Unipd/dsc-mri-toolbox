function [fwhm_map]=DSC_mri_fwhm(conc,mask,options)
conc_vect = vol2mat(conc,mask);
fwhm_vect = zeros(size(conc_vect,1),1);

for r=1:size(conc_vect,1)
    
    try
        temp = fwhm(options.time,conc_vect(r,:));
    catch
        temp = 0;
    end
    fwhm_vect(r) = temp;
end


fwhm_map = zeros(size(mask));
fwhm_map(logical(mask)) = fwhm_vect;
end

