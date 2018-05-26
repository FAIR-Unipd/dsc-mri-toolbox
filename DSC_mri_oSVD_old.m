function [res_osvd]=DSC_mri_oSVD(conc,aif,mask,options)
% ultima modifica: Marco Castellaro 10/03/2013

%Funzione del pacchetto DSC_mri - DSC_mri_cbf
%Autore: Castellaro Marco - Università di Padova - DEI
%
%Calcola le mappe parametriche di Cerebral Blood Flow (CBF) per un soggetto
%e per la deconvoluzione utilizza il metodo della SINGULAR VALUE
%DECOMPOSITION versione BLOCK-CIRCULANT con Indice di oscillazione
%
%Parametri in ingresso - conc (Matrice 4D), contiene gli andamenti delle
%                        concentrazioni DSC di tutti i voxel.
%                      - aif, andamento delle concentrazioni nel sito
%                        scelto come arteriale
%                      - mask (Matrice 3D), contiene la matrice utilizzata
%                        per mascherare il volume cerebrale da analizzare
%Options è la sruct che contiene le opzioni del metodo, quelli
%significativi sono:
%
%options.deconv.osvd.OI -       Obscillation index - default 0.035
%
%options.deconv.SVD.residual - se a 1 produce in uscita anche la matrice 4D
%                              dei residui, altrimenti non vengono
%                              calcolati

if options.display > 1
    disp('    Method: oSVD');
    disp(['    OI: ' num2str(options.deconv.oSVD.OIthres,'%1.3f') ]);
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
OIthres=options.deconv.oSVD.OIthres;    % threshold del 10% con in Ostergaard e Calamante
OIcounter=options.deconv.oSVD.OIcounter;

% 3) Applico la Ginv per calcolare la funzione residuo e il CBF di ciascun
%    voxel.
res_osvd.map=zeros(options.nR,options.nC,options.nS);
if options.deconv.oSVD.residual
    res_osvd.residual=zeros(options.nR,options.nC,options.nS,nTpad);
end

[U,S,V]=svd(G);
W=diag(1./diag(S));

if options.waitbar
    hw_osvd=waitbar(0,'Computing CBF by oSVD');
end
for r=1:options.nR
    if options.waitbar
        waitbar((r-1)/(options.nR),hw_osvd);
    end
    for c=1:options.nC
        for s=1:options.nS
            if mask(r,c,s)
                % Calcolo la funzione residuo
                vettConc=zeros(nTpad,1);
                vettConc(1:options.nT)=reshape(conc(r,c,s,:),options.nT,1);
                
                for threshold=5:5:95;
                    W(S<threshold/100*S(1,1))=0;
                    vettRes=V*W*(U'*vettConc);
                    O=0;L=length(vettRes);
                    for j1=3:L
                        O=O+abs(vettRes(j1)-2*vettRes(j1-1)+vettRes(j1-2));
                    end
                    OI=1/L*1/max(vettRes)*O;
                    
                    if OI<OIthres(OIcounter)
                        break
                    end
                end
                
                res_osvd.map(r,c,s)=max(vettRes);
                if options.deconv.oSVD.residual
                    res_osvd.residual(r,c,s,:)=vettRes;
                end
            end
        end
    end
end
if options.waitbar
    delete(hw_osvd)
end


end