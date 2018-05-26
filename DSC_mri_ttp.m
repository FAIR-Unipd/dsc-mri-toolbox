function [ttp]=DSC_mri_ttp(conc,mask,options)
[~,ttp]=max(conc,[],4);
ttp(not(logical(mask))) = 0;


end

