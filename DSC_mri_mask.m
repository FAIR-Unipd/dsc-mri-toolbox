function [mask]=DSC_mri_mask(volumes,options)
% ultima modifica: Marco Castellaro 08/07/2010

%Funzione del pacchetto DSC_mri - DSC_mri_mask
%Autore: Castellaro Marco - Università di Padova - DEI
%
%Calcola le maschere per esami DSC-MRI.
%
%Parametri in ingresso: volumes (Matrice 4D) che contiene gli andamenti del
%segnale DSC di tutti i voxel.
%Options è la sruct che contiene le opzioni del metodo, quelli
%significativi sono:
%
%options.mask.npixel: rappresenta il numero di pixel minimi di una
%                     componente connessa che è utilizzata come soglia
%                     per escludere dall'immagine lo scalpo e le zone
%                     adiacenti all'esterno dell'encefalo
%
%options.display - livello 1 Mostra l'avanzamento dell'elaborazione
%                - livello 2 Mostra le maschere e informazioni sulla soglia
%                  e sulle intensità delle immagini da mascherare
%
%Parametri in uscita: struttura mask, che contiene
% - aif: maschera ottimizzata per la ricerca della funzione di ingresso arteriale
% - data: maschera ottimizzata per il mascheramento dell'intero encefalo
% - threshold: soglia calcolata e fornita in uscita

if options.display > 0
    disp(' ')
    disp('Masking data... ');
end

nbin=100;

volume_sum=sum(volumes,4);
mask.data=zeros(size(volume_sum));
mask.aif=zeros(size(volume_sum));

[prob,intensity]=hist(volume_sum(1:options.nR*options.nC*options.nS),nbin);

[~,ind_max]=max(prob);
if (prob(ind_max+1)/prob(ind_max)) > 0.3 
    ind_max=ind_max-1;
end

temp=prob(ind_max+1:end);
clear prob; prob=temp;

temp=intensity(ind_max+1:end);
clear intensity; intensity=temp;

%Metodo Marco v 2.0 - fit 2 gaussiane

f = fittype('gauss2');
g1=inline('a1.*exp(-((x-b1)./c1).^2)');

opt_fit = fitoptions('gauss2');
gfit = fit(double(intensity)',double(prob)',f,opt_fit);

[mask.threshold,~]=curveintersect(intensity,g1(gfit.a1,gfit.b1,gfit.c1,intensity), ... 
                                  intensity,g1(gfit.a2,gfit.b2,gfit.c2,intensity));


% % Metodo Denis per identificare la soglia
% [maxValue,maxPos1]=max(prob); % Individua il primo picco delle intensità
% exitFlag=true;
% maxPos2=maxPos1+1;
% while exitFlag % Cerca il secondo picco delle intensità
%     [maxValue,maxPos]=max(prob(maxPos2:nbin));
%     if (maxPos~=1)||(maxPos2==nbin)
%         exitFlag=false;
%         maxPos2=maxPos2+maxPos-1;        
%     else
%         maxPos2=maxPos2+1;
%     end
% end
% 
% [minValue,minPos]=min(prob(maxPos1:maxPos2)); % Individua il minimo tra i due picchi di intensità
% mask.threshold=intensity(minPos+maxPos1-1);



% % Metodo Marco per identificare la soglia
% t0=0.1*nbin+1;
% 
% der_prob=(prob(2:end)-prob(1:end-1))./(prob(2)-prob(1));
% der_probS = smooth(der_prob,t0,'moving');
% der2_prob=(der_probS(2:end).*der_probS(1:end-1)) <0;
% 
% der2_prob=find(der2_prob>0);
% 
% mask.threshold=intensity(der2_prob(1)+2);
% 
% 
if options.display > 0
    if size(mask.threshold,1) > 1
        mask.threshold=mask.threshold(1);
        disp(['   Threshold (multiple intercept - 1) : ' num2str(mask.threshold)]);
    else
        disp(['   Threshold: ' num2str(mask.threshold)]);
    end
    if options.display > 1
        
        hf_hist=figure();
        hist(volume_sum(1:options.nR*options.nC*options.nS),nbin)
        hold on
        line([mask.threshold mask.threshold],[0 1.05*max(prob)], ...
            'Color','g','LineStyle','--');
        plot(intensity,gfit(intensity),'r')
        ylim([0 1.05*max(prob)])
    end
end

mask.aif=volume_sum>mask.threshold;
hf_mask=zeros(options.nS,1);

for s=1:options.nS
    %copro eventuali "buchi" creati dalla sogliatura
    temp = mask.aif(:,:,s);
    
    %elimino le componenti connesse minori e lascio intatte quelle maggiori
    %di #options.mask.pixel
    
    CC = bwconncomp(temp);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    
    for j=1:size(numPixels,2)
        
        if numPixels(idx) > numPixels(j) && numPixels(j) < options.mask.npixel
            
            temp(CC.PixelIdxList{j}) = 0;
        end
    end
    
    mask.data(:,:,s)=imfill(temp,'holes');
    
    if (options.display > 2)||((options.display >1) && (s == round(0.5*options.nS)))
        hf_mask(s)=figure();
        
        subplot(121)
        imagesc(volume_sum(:,:,s))
        colormap('gray')
        title(['Slice ' num2str(s) '/' num2str(options.nS) ' - Masked Data for AIF selection'])
        
        B = bwboundaries(temp,4);
        hold on
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
        end
        
        subplot(122)
        imagesc(volume_sum(:,:,s))
        
        B = bwboundaries(mask.data(:,:,s),4);
        hold on
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
        end
        colormap('gray')
        title(['Slice ' num2str(s) '/' num2str(options.nS) ' - Masked Data'])
        
    end
end
mask.aif=mask.aif.*mask.data;

if options.display > 0
    if options.display > 2
        pause
        for cont_h=1:length(hf_mask)
            try
                close(hf_mask(cont_h))
            end          
        end
        for cont_h=1:length(hf_hist)
            try
                close(hf_hist(cont_h))
            end
            
        end
    end
end

mask.gfit=gfit;
mask.intensity=intensity;
end

