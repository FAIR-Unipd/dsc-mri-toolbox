function [res_svd]=DSC_mri_SVD(conc,aif,mask,options)
%Funzione del pacchetto DSC_mri - DSC_mri_cbf
%Autore: Castellaro Marco - Universit� di Padova - DEI
%
%Calcola le mappe parametriche di Cerebral Blood Flow (CBF) per un soggetto
%e per la deconvoluzione utilizza il metodo della SINGULAR VALUE
%DECOMPOSITION con troncamento
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
%                               all'autovalore massimo, in Ostergaard e
%                               Calamante et al. � fissato al 20%
%
%options.deconv.SVD.residual - se a 1 produce in uscita anche la matrice 4D
%                              dei residui, altrimenti non vengono
%                              calcolati

if options.display > 1
    disp('    Method: SVD');
    disp(['    Threshold: ' num2str(100*options.deconv.SVD.threshold,'%3d') '%']);
end

% 1) Creo la matrice G
aifVett=zeros(options.nT,1);
aifVett(1)=aif(1);
aifVett(options.nT)=aif(options.nT);

for k=2:(options.nT-1)
    aifVett(k)=(aif(k-1)+4*aif(k)+aif(k+1))/6;
end

G=toeplitz(aifVett,[aifVett(1) zeros(1,options.nT-1)]);

% 2) Applico la SVD per calcolare la G inversa
[U,S,V]=svd(G);

eigenV=diag(S);
threshold=options.deconv.SVD.threshold*max(eigenV);    % threshold del 10% in Ostergaard e Calamante
newEigen=zeros(size(eigenV));
for k=1:length(eigenV);
    if eigenV(k)>=threshold;
        newEigen(k)=1/eigenV(k);
    end
end

Ginv=V*diag(newEigen)*(U');

% 3) Applico la Ginv per calcolare la funzione residuo e il CBF di ciascun
%    voxel.
res_svd.map=zeros(options.nR,options.nC,options.nS);
if options.deconv.SVD.residual
	res_svd.residual=zeros(options.nR,options.nC,options.nS,options.nT);
end

if options.waitbar
    hw_svd=waitbar(0,'Computing CBF by SVD');
end
for r=1:options.nR
    if options.waitbar
        waitbar((r-1)/(options.nR),hw_svd);
    end
    for c=1:options.nC
        for s=1:options.nS
            if mask(r,c,s)
                % Calcolo la funzione residuo
                
                vettConc=reshape(conc(r,c,s,:),options.nT,1);
                vettRes=(1/options.tr)*Ginv*vettConc;
                
                res_svd.map(r,c,s)=max(abs(vettRes));
                if options.deconv.SVD.residual
                    res_svd.residual(r,c,s,:)=vettRes;
                end
                
            end
            
        end
    end
end

if options.waitbar
    delete(hw_svd)
end
end