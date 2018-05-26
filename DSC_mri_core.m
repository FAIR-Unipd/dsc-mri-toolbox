function [cbv,cbf,mtt,ttp,mask,aif,conc,s0,cbv_last,fwhm,K1_map,K2_map,K1_CV_map,K2_CV_map]=DSC_mri_core(volumes,te,tr,optionsIN,aif,mask)
% ultima modifica: Marco Castellaro 26/05/2018

addpath('utils')

%----- CHECK OF INPUT PARAMETERS - PART 1 ----------------------
% Check that input and output has the minimum number of mandatory parameters 
% Input
if nargin<1
    error('Input data required.')
end

% Output
if nargout < 3
    error('Required almost tree exiting parameters: cbf,cbv and mtt.')
end

% Set Default parameters in case of not providing TE and TR
switch nargin
    case 1
        tr=1;
        te=1;
    case 2
        tr=1;
end

% Check if options are provided, otherwise load default parameters
if nargin < 4
    % Load default options
    options=DSC_mri_getOptions;  
else
    % Use provided options
    options=optionsIN;
    clear optionsIN;
end

%----- INTRO --------------------------------------------------------------
if options.display>0
    disp(' ');
    disp(' ');
    disp(' ');
    disp('  _____   ____    _____         _    _   ____    ___ ')
    disp(' |  _  \ /  __|  /  ___|       | \  / | |  _  \ |   |')
    disp(' | | | | | |__   | |      ___  |  \/  | | |_| |  | | ')
    disp(' | | | | \___ \  | |     |___| | |\/| | |     /  | | ')
    disp(' | |_| |  ___| | | |___        | |  | | | |\ \   | | ')
    disp(' |_____/ |____/  \_____|       |_|  |_| |_| \_\ |___|')
    disp(' ')
    disp('                       TOOLBOX                       ')
    disp('          by Denis Peruzzo & Marco Castellaro')
    disp(' ')
end

%----- INSERIMENTO DIMENSIONI VETT. OPZIONI ------------------------------
if options.display > 0
    disp('Checking data...');
end
[nR,nC,nS,nT]=size(volumes);
options.nR=nR;
options.nC=nC;
options.nS=nS;
options.nT=nT;

%----- INSERIMENTO TEMPI NEL VETT. OPZIONI -------------------------------
options.te=te;
options.tr=tr;
options.time=0:tr:(nT-1)*tr;

if options.display > 0
    disp('                      DATA SIZE ');
    disp(['                     Rows - ' num2str(options.nR)]);
    disp(['                  Columns - ' num2str(options.nC)]);
    disp(['                   Slices - ' num2str(options.nS)]);
    disp(['                  Samples - ' num2str(options.nT)]);
    disp(['                Echo time - ' num2str(options.te) ' s']);
    disp(['          Repetition time - ' num2str(options.tr) ' s']);
end

%----- CHECK OF INPUT PARAMETERS - PART 2 ----------------------
% Check wether the AIF is column of row vector
if nargin > 4
    if not(isempty(aif))
        if not(nT==size(aif.conc,1))
            if not(nT==size(aif.conc,2))
                if size(aif.conc,1)==0 && size(aif.conc,2)==0
                    options.aif=1;
                else
                    error('AIF dimensions are not compatible to samples data, please provide a compatible AIF.')
                end
            end
        end
    end
end


%----- MASK CALCULATION --------------------------------------------
if nargin < 5
    if not(options.conc)
        [mask]=DSC_mri_mask(volumes,options);
    else
        error('In case of concentration as first parameter mask is required. Use DSC_mri_core(conc, ... ,MASK).')
    end
elseif size(mask,1)==0
    if not(options.conc)
        [mask]=DSC_mri_mask(volumes,options);
    end
else
    if options.aif.enable
        [mask]=DSC_mri_mask_only_aif(volumes,mask,options);
        % Check mask dimension 
        if not(nR==size(mask.data,1))
            error('Rows mask not equal to data.')
        elseif not(nC==size(mask.data,2))
            error('Columns mask not equal to data.')
        elseif not(nS==size(mask.data,3))
            error('Slices mask not equal to data.')
       
        end
    else
    
        % Check mask dimension 
        if not(nR==size(mask,1))
            error('Rows mask not equal to data.')
        elseif not(nC==size(mask,2))
            error('Columns mask not equal to data.')
        elseif not(nS==size(mask,3))
            error('Slices mask not equal to data.')
        else
            temp=mask;
            clear mask;
            mask.data=temp;
        end
    end
end

%----- CALCOLO DELLE CONCENTRAZIONI E DI S0 ------------------------------
if options.conc
    conc=volumes;
    options.aif.nSlice=aif.sliceAIF;
    %----- CALCOLO SI s0 - OPZIONALE -----------------------------------------
elseif nargout > 7
    [conc,s0,bolus]=DSC_mri_conc(volumes,mask.data,options);
    conc = real(conc);
else
    [conc,s0,bolus]=DSC_mri_conc(volumes,mask.data,options);
    conc = real(conc);
end

%----- ESTRAZIONE DELL'AIF -----------------------------------------------
if options.aif.enable
    
    [aif]=DSC_mri_aif(conc,mask.aif,options);
end


%----- CALCOLO DELLE MAPPE -----------------------------------------------
if options.display>0
    disp(' ')
    disp('Calculating perfusion maps...')
end

    
%----- CALCOLO DEL CBV ---------------------------------------------------
[cbv]=DSC_mri_cbv(conc,aif.fit.gv,mask,options);

%----- CALCOLO DEL CBV (leackage monitoring) ----
[cbv_last,K1_map,K2_map,K1_CV_map,K2_CV_map]=DSC_mri_cbv_lc(conc,aif.fit.gv,mask,bolus,options);

%----- CALCOLO DEL CBF ---------------------------------------------------
[cbf]=DSC_mri_cbf(conc,aif.fit.gv,mask,options);

%----- CALCOLO DEL MTT ---------------------------------------------------
[mtt]=DSC_mri_mtt(cbv,cbf,options);

%----- CALCOLO DEL TTP - OPZIONALE ---------------------------------------
if nargout > 3
    [ttp]=DSC_mri_ttp(conc,mask.data,options);
    [fwhm]=DSC_mri_fwhm(conc,mask.data,options);
end


%----- PREPARAZIONE DELLE USCITE -----------------------------------------
%Viene fornita in uscita solo la matrice con le maschere da applicare sulle
%mappe parametriche
mask=mask.data;

end