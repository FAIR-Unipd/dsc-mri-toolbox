function [cbv_corrected,K1_map,K2_map,K1_CV_map,K2_CV_map]=DSC_mri_cbv_lc(conc,aif,mask,bolus,options)
%Funzione del pacchetto DSC_mri - DSC_mri_cbv
%Autore: Castellaro Marco - Universit� di Padova - DEI
%
%Calcola le mappe parametriche di Cerebral Blood Volume corrette per leakage per un soggetto
%
%Parametri in ingresso: conc (Matrice 4D) che contiene gli andamenti delle
%concentrazioni DSC di tutti i voxel.
%Options � la sruct che contiene le opzioni del metodo, quelli
%significativi sono:
%
%options.time - Rappresenta il vettore dei tempi dell'esame DSC (ogni
%               campione rappresenta l'acquisizione di un intero volume cerebrale
%
%options.par.kh - Parametro che rappresenta la dipendenza dall'ematocrito
%                 del Cerebral Blood Volume (CBV), per default � settato ad
%                 uno, in questo caso si ottengono stime RELATIVE del
%                 parametro CBV.
%
%options.par.rho - Parametro che rappresenta la dipendenza dalla densit�
%                  del sangue del Cerebral Blood Volume (CBV), per default
%                  � settato ad uno, in questo caso si ottengono stime
%                  RELATIVE del parametro CBV.
%
%Parametri in uscita:
%cbv - (matrice 3D) che contiene al suo interno la mappa parametrica
%      calcolata


if options.display > 0
    disp('   CBV');
end

cbv=zeros(options.nR,options.nC,options.nS);

if options.waitbar
    hw_cbv=waitbar(0,'Calculating CBV last');
end

% We next estimated R2*(t) by averaging R2*(t) for all pixels within the mask that did not demonstrate signal intensity enhancement (averaged over the final 10 time points) greater than 1 SD above that pixel’s
% average baseline and estimated by using trapezoidal integration over the 120 acquired time points
% 


SD = vol2mat (conc,mask.data);
SD = nanstd(SD(:,1:floor(mean(bolus(:)))),[],2);
SD_map = zeros(size(cbv));
SD_map(mask.data) = SD;

AVG = vol2mat (conc,mask.data);
AVG = nanmedian(AVG(:,end-10:end),2);
AVG_map = zeros(size(cbv));
AVG_map(mask.data) = AVG;


mask_not_enhancing = and(mask.data, not(abs(AVG_map) > 2.*SD_map));
conc(isinf(conc)) = 0;
R2star_AVG_not_enhancing = nanmean(vol2mat(conc,mask_not_enhancing))';


Delta_R2star = vol2mat(conc,mask.data);

phat = zeros(size(Delta_R2star,1),2);
CVp = zeros(size(Delta_R2star,1),2);
A = [-cumtrapz(options.time,R2star_AVG_not_enhancing)  R2star_AVG_not_enhancing];
sigmaphat = inv(A'*A);

bolus_min = min(bolus(:));

for v=1:size(Delta_R2star,1)
    Delta_R2star_vett = Delta_R2star(v,:)';

    sigma2 = (std(Delta_R2star_vett(1:bolus_min))).^2;
    
    B = Delta_R2star_vett;
    phat(v,:) = A\B;
    temp = sigma2.*sigmaphat;
    CVp(v,:) = 100.*sqrt([temp(1,1) temp(2,2)])./abs(phat(v,:));
end

K2_vett  = phat(:,1);
K1_vett  = phat(:,2);

K2_CV_vett  = CVp(:,1);
K1_CV_vett  = CVp(:,2);

K2_map  = zeros(size(cbv));
K1_map  = zeros(size(cbv));

K2_map(mask.data) = K2_vett;
K1_map(mask.data) = K1_vett;

K2_CV_map = zeros(size(cbv));
K1_CV_map = zeros(size(cbv));


K2_CV_map(mask.data) = K2_CV_vett;
K1_CV_map(mask.data) = K1_CV_vett;

for s=1:options.nS
    cbv(:,:,s)=mask.data(:,:,s).* ...
        (trapz(options.time,conc(:,:,s,:),4));
    if options.waitbar
        waitbar(s/options.nS,hw_cbv)
    end
end

cbv_corrected = (cbv + (K2_CV_map<100).*K2_map.*trapz(options.time,cumtrapz(options.time,R2star_AVG_not_enhancing)))./trapz(options.time,aif);

% figure(); 
% subplot(131); 
% imagesc(imrotate(cbv(:,:,14),90),[0 1]);
% colorbar
% subplot(132); 
% imagesc(imrotate(cbv_corrected(:,:,14),90),[0 10]);
% colorbar
% subplot(133); 
% imagesc(imrotate(K2_map(:,:,14),90),[-0.4 0.4]);
% colorbar
% colormap('jet')



if options.waitbar
    delete(hw_cbv);
end
end
