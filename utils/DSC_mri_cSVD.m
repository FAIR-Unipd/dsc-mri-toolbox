function [res_csvd]=DSC_mri_cSVD(conc,aif,mask,options)
% ultima modifica: Denis Peruzzo 07/06/2010

%Funzione del pacchetto DSC_mri - DSC_mri_cbf
%Autore: Castellaro Marco - Universit� di Padova - DEI
%
%Calcola le mappe parametriche di Cerebral Blood Flow (CBF) per un soggetto
%e per la deconvoluzione utilizza il metodo della SINGULAR VALUE
%DECOMPOSITION versione BLOCK-CIRCULANT con troncamento
%
%Parametri in ingresso - conc (Matrice 4D), contiene gli andamenti delle
%                        concentrazioni DSC di tutti i voxel.
%                      - aif, andamento delle concentrazioni nel sito
%                        scelto come arteriale
%                      - mask (Matrice 3D), contiene la matrice utilizzata
%                        per mascherare il volume cerebrale da analizzare
%Options � la sruct che contiene le opzioni del metodo, quelli
%significativi sono:
%
%options.deconv.svd.threshold - soglia di troncamento percentuale riferita
%                               all'autovalore massimo, in Wu et al.
%                               Calamante � fissato al 10%
%
%options.deconv.SVD.residual - se a 1 produce in uscita anche la matrice 4D
%                              dei residui, altrimenti non vengono
%                              calcolati

if options.display > 1
    disp('    Method: cSVD');
    disp(['    Threshold: ' num2str(100*options.deconv.cSVD.threshold,'%3d') '%']);
end

% 1) Creo la matrice G
nTpad=2*options.nT;
columnG=zeros(nTpad,1);
columnG(1)=aif(1);
columnG(options.nT)=(aif(options.nT-1)+4*aif(options.nT))/6;
columnG(options.nT+1)=aif(options.nT)/6;
for k=2:(options.nT-1)
    columnG(k)=(aif(k-1)+4*aif(k)+aif(k+1))/6;
end
rowG=zeros(1,nTpad);
rowG(1)=columnG(1);
for k=2:nTpad
    rowG(k)=columnG(nTpad+2-k);
end

G=toeplitz(columnG,rowG);

% 2) Applico la SVD per calcolare la G inversa
[U,S,V]=svd(G);

eigenV=diag(S);
threshold=options.deconv.cSVD.threshold*max(eigenV);    % threshold del 10% con in Ostergaard e Calamante
newEigen=zeros(size(eigenV));
for k=1:length(eigenV);
    if eigenV(k)>=threshold;
        newEigen(k)=1/eigenV(k);
    end
end

Ginv=V*diag(newEigen)*(U');

% 3) Applico la Ginv per calcolare la funzione residuo e il CBF di ciascun
%    voxel.
res_csvd.map=zeros(options.nR,options.nC,options.nS);
if options.deconv.SVD.residual
    res_csvd.residual=zeros(options.nR,options.nC,options.nS,nTpad);
end

if options.waitbar
    hw_csvd=waitbar(0,'Computing CBF by cSVD');
end
for r=1:options.nR
    if options.waitbar
        waitbar((r-1)/(options.nR),hw_csvd);
    end
    for c=1:options.nC
        for s=1:options.nS
            if mask(r,c,s)
                % Calcolo la funzione residuo
                vettConc=zeros(nTpad,1);
                vettConc(1:options.nT)=reshape(conc(r,c,s,:),options.nT,1);
                vettRes=(1/options.tr)*Ginv*vettConc;
                
                res_csvd.map(r,c,s)=max(abs(vettRes));
                if options.deconv.cSVD.residual
                    res_csvd.residual(r,c,s,:)=vettRes;
                end
            end
        end
    end
end
if options.waitbar
    delete(hw_csvd)
end


end