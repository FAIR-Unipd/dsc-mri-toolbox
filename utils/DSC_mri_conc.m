function [conc,S0map,bolus]=DSC_mri_conc(volumes,mask,options)
% ultima modifica: Denis Peruzzo 07/06/2010

%Funzione del pacchetto DSC_mri - DSC_mri_conc
%Autore: Castellaro Marco - Universit� di Padova - DEI
%Data:
%
%Calcola le concentrazioni e le mappe di s0 in esami DSC-MRI.
%
%Parametri in ingresso:
% - volumes (matrice 4D) che contiene gli andamenti del
%   segnale DSC di tutti i voxel.
% - mask (matrice 3D) contiene la matrice per mascherare i voxel non di
%   interesse per lo studio
%
%Options � la sruct che contiene le opzioni del metodo, quelli
%significativi sono:
%
% par.kvoi - Costante di proporzionalit� per il calcolo della
%            concentrazione di tracciante nel VOI, per default �
%            considerata non nota e posta a 1.
%
% S0 - serie di parametri per individuare la soglia di calcolo di S0
% S0.nSamplesMin  - n� di campioni che considero sicuramente acquisiti
%                   prima dell'iniezione
% S0.nSamplesMax  - n� di campioni dopo il quale mi fermo in ogni caso
% S0.thresh;      - Aggiungo un campione se la sua diff. dalla media �
%                   minore della soglia
%
%options.display - livello 1 Mostra l'avanzamento dell'elaborazione
%                - livello 2 Mostra le mappe di S0 e il segnale medio su
%                  cui lo si � stimato
%
%Parametri in uscita:
% - conc: matrice 4D delle concentrazioni
% - S0map: matrice 3D degli S0


if options.display > 0
    disp(' ')
    disp('Calculating concentration...');
end

[S0map,bolus]=DSC_mri_S0(volumes,mask,options);

conc=zeros(size(volumes));

ind=find(mask);
k=options.nR*options.nC*options.nS;
if options.waitbar
    hw=waitbar(0,'Calculating concentration...');
end

for t=1:options.nT
    if options.waitbar
        waitbar((t-1)/options.nT,hw)
    end
    step1=volumes(ind+k*(t-1))./S0map(ind);
    %step2=step1.*(step1<1);
    %step3=step2+(step2==0);
    conc(ind+k*(t-1))=-(options.par.kvoi/options.te).*log(step1);
end

if options.waitbar
    delete(hw);
end

end

%%
function [S0map,bolus]=DSC_mri_S0(volumes,mask,options)
%Funzione del pacchetto DSC_mri - DSC_mri_S0
%Autore: Castellaro Marco - Universit� di Padova - DEI

% La funzione calcola l'istante di iniezione del bolo e l'S0 a partire dai
% dati.
% 1) sull'andamento medio calcolo l'istante di iniezione del bolo: calcolo
%    la media dei primi n campioni e aggiungo l'n+1 se la sua differenza
%    percentuale dalla media � inferiore ad una data soglia.
% 2) calcolo S0 come la media dei primi n campioni per tutti i voxel.

% Definizione dei parametri
nSamplesMin=options.S0.nSamplesMin;
% n� di campioni che considero sicuramente acquisiti prima dell'iniezione.
nSamplesMax=options.S0.nSamplesMax;
% n� di campioni dopo il quale mi fermo in ogni caso.
thresh=options.S0.thresh;
% Aggiungo un campione se la sua diff dalla media � minore della soglia.

% 1) calcolo dell'istante di iniezione del bolo
% 1.1) calcolo dell'andamento medio

mean_signal=zeros(size(options.time));

if options.waitbar
    hw=waitbar(0,'Retriving S0...');
end
for s=1:options.nS
    for t=1:options.nT
        if options.waitbar
            waitbar((options.nT*(s-1)+(t-1))/(options.nT*options.nS),hw);
        end
        indMask=find(mask(:,:,s));
        mean_signal(s,t)=mean(volumes(indMask+(options.nR*options.nC*options.nS)*(t-1)+(options.nS*(s-1))));
        
        clear('temp','indMask');
    end
end

for s=1:options.nS
    % 1.2) calcolo dell'istante di iniezione del bolo
    ciclo=true;
    pos=nSamplesMin;
    while ciclo
        mean_val=mean(mean_signal(s,1:pos));
        if abs((mean_val-mean_signal(s,pos+1))/mean_val)<thresh
            pos=pos+1;
        else
            ciclo=false;
            pos=pos-1; % Scelta conservativa, non considero
            % l'ultimo campione prima dell'iniezione
        end
        if pos==nSamplesMax
            ciclo=false;
            pos=pos-1;
        end
    end
    
    if (options.display > 2)||((s==round(0.5*options.nS))&&options.display > 1)
        hf_s0(s)=figure();
        figure(hf_s0(s))
        subplot(1,2,1)
        plot(options.time,mean_signal(s,:),'b+-')
        hold on
        plot(options.time(pos),mean_signal(s,pos),'ro')
        %legend('Mean Signal','Final sample time used to calc S0')
        xlabel(['S0 computed from first ' num2str(pos) ' samples.'])
        xlim([options.time(1) options.time(end)])
        title(['Slice ' num2str(s) '/' num2str(options.nS) ' - bolus arrival time'])
    end
    
    % 2) Calcolo di S0
    S0map(:,:,s)=mask(:,:,s).*mean(volumes(:,:,s,1:pos),4);
   
    
    if (options.display > 2)||((s==round(0.5*options.nS))&&options.display > 1)
        figure(hf_s0(s))
        subplot(1,2,2)
        imagesc(S0map(:,:,s)),colorbar
        title('S0')
    end
    bolus(s) = pos;
end
if options.display > 2
    pause()
    for cont_h=1:length(hf_s0)
        try
            close(hf_s0(cont_h))
        end
    end
end
if options.waitbar
    try
        delete(hw);
    end
end
end

