function [aif]=DSC_mri_aif(conc,mask,options)
% ultima modifica: Denis Peruzzo 08/06/2010

%Funzione del pacchetto DSC_mri - DSC_mri_aif
%Autore: Denis Peruzzo - Universit� di Padova - DEI
%
%Individua la AIf per esami DSC-MRI. Il metodo � pensato per individuare la
%middle cerebral artery (MCA) nelle slice subito superiori al corpo
%calloso.
%
%Parametri in ingresso: 
% - conc (Matrice 4D) che contiene gli andamenti di concentrazione di tutti
%   i voxel
% - mask (matrice 3D) contiene la maschera di ogni slice. NB: la maschera 
%   non � quella utilizzata per il calcolo delle concentrazioni, ma una sua
%   versione restrittiva (non � stato utilizzata la funzione fill)
% - options � la sruct che contiene le opzioni del metodo, quelli 
%   significativi sono:
%
% options.aif.enable: flag che abilita la ricerca dell'AIF (default=1)
%
% options.aif.semiasseMaggiore: permette di individuare l'area di ricerca
%                               della AIF. La zona di ricerca elittica
%                               presenta il seiasse maggiore pari a due
%                               volte la porzione indicata da questa
%                               opzione della dimensione dell'encefalo
%                               (default=0.35)
% 
% options.aif.semiasseMinore: come il semiasse maggiore, ma relativo 
%                             all'altro semiasse della regione di ricerca
%                             (default=0.15)
%
% options.aif.pArea: percentuale di voxel candidati esclusi sulla base
%                    dell'area sotto la curva (default=0.4)
% 
% options.aif.pTTP: percentuale di voxel candidati esclusi sulla base del
%                   time to peak (default=0.4)
% 
% options.aif.pReg: percentuale di voxel candidati esclusi sulla base
%                   dell'irregolarit� della curva (default=0.05)
% 
% options.aif.diffPicco: 0.0400
% 
% options.aif.nVoxelMax: numero massimo di voxel arteriali accettati
%                        (default=6)
% 
% options.aif.nVoxelMin: numero minimo di voxel arteriali accettati
%                        (default=4)
%
%Parametri in uscita: struttura aif, che contiene
% - ROI

%%
if options.display > 0
    disp(' ')
    disp('AIF extraction...');
end

%Selezione slice AIF
if ( options.aif.nSlice <= 0 )||( options.aif.nSlice > options.nS )
    options.aif.nSlice=DSC_mri_aif_slice_selection_figure(conc,options);
end


[aif_old]=estraiAIF(reshape(conc(:,:,options.aif.nSlice,:),options.nR,options.nC,options.nT),mask(:,:,options.aif.nSlice),options);
aif=aif_old;
aif.sliceAIF=options.aif.nSlice;

if options.qr.enable
    
    aif.conc=(options.qr.a*aif_old.conc+options.qr.b*(aif_old.conc.^2))/options.qr.r;
    if options.display > 1
        hf.aif=figure();
        subplot(121)
        plot(options.time,aif_old.conc)
        title('AIF whit linear \Delta R_2^* relationship.')
        subplot(122)
        plot(options.time,aif.conc)
        title('AIF whit quadratic \Delta R_2^* relationship.')
        pause()
        close(hf.aif)
        
    end
end




%%
function [AIF]=estraiAIF(AIFslice,mask,options)
% Estrae la AIF dalla slice fornita in ingresso.
% 1) Individua la regione contenente la AIF
% 2) Decimazione dei voxel candidati
% 3) Applica l'algoritmo di cluster gerarchico per individuare i voxel
%    arteriali
% 4) Prepara l'output

% Preparazione variabili accessorie e dei parametri
semiasseMag = options.aif.semiasseMaggiore;
semiasseMin = options.aif.semiasseMinore;
pArea       = options.aif.pArea;
pTTP        = options.aif.pTTP;
pReg        = options.aif.pReg;
nVoxelMax   = options.aif.nVoxelMax;
nVoxelMin   = options.aif.nVoxelMin;
diffPicco   = options.aif.diffPicco;


[nR,nC,nT]=size(AIFslice);
maschera=mask;    % MARCO

% Preparo l'immagine per le eventuali visualizzazioni
immagine.img=sum(AIFslice,3);
vettImmagine=sort(immagine.img(1:nR*nC));
immagine.bound=[0 vettImmagine(round(0.95*nR*nC))];
clear vettImmagine


% 1) INDIVIDUAZIONE DELLA ROI CONTENENTE LA AIF
% 1.1) Individuazione degli estremi della maschera
if options.display > 0
    disp('   Brain bound detection')
end
ciclo=true;
r=1;
while ciclo
    if sum(maschera(r,:))~=0
        minR=r;
        ciclo=false;
    else
        r=r+1;
    end
end
ciclo=true;
r=nR;
while ciclo
    if sum(maschera(r,:))~=0
        maxR=r;
        ciclo=false;
    else
        r=r-1;
    end
end

ciclo=true;
c=1;
while ciclo
    if sum(maschera(:,c))~=0
        minC=c;
        ciclo=false;
    else
        c=c+1;
    end
end
ciclo=true;
c=nC;
while ciclo
    if sum(maschera(:,c))~=0
        maxC=c;
        ciclo=false;
    else
        c=c-1;
    end
end

if options.display > 2
    hf.mask=figure();
    imagesc(maschera)
    hold on
    plot([1 options.nC],[minR minR],'g-',[1 options.nC],[maxR maxR],'g-',[minC minC],[1 options.nR],'g-',[maxC maxC],[1 options.nR],'g-')
    xlabel(['Brain bound: rows (' num2str(minR) '-' num2str(maxR) ') - columns(' num2str(minC) '-' num2str(maxC) ')'])
    title('AIF extraction - mask and bounds')
    set(gca,'xtick',[],'ytick',[])
    axis square
end

% 1.2) Disegno della ROI
if options.display > 0
    disp('   Definition of the AIF extration searching area')
end
centro(2)=0.5*(minR+maxR); % Coordinata Y del centro (calcolata sulle righe)
centro(1)=0.5*(minC+maxC); % Coordinata X del centro (calcolata sulle colonne)

% Il semiasse maggiore � lungo la direzione antero-posteriore e quindi da
% sinistra a destra delle immagini.
semiasseB=semiasseMag.*(maxC-minC);
semiasseA=semiasseMin.*(maxR-minR);
ROI=zeros(nR,nC); % maschera contenente i voxel della ROI
for r=1:nR
    for c=1:nC
        if ((r-centro(2))^2/(semiasseA^2)+(c-centro(1))^2/(semiasseB^2))<=1
            ROI(r,c)=1;
        end
    end
end
ROI=ROI.*maschera; % Dei voxel della ROI mantengo solo quelli presenti anche nella maschera
ROIiniziale=ROI;

xROI=centro(1)-semiasseB:0.01:centro+semiasseB;
nL=length(xROI);
xROI(2*nL)=0;
yROI=zeros(1,2*nL);
for k=1:nL
    yROI(k)=semiasseA*((1-((xROI(k)-centro(1))^2)/(semiasseB^2))^0.5)+centro(2);
end
for k=1:nL
    xROI(nL+k)=xROI(nL-k+1);
    yROI(nL+k)=-semiasseA*((1-((xROI(nL+k)-centro(1))^2)/(semiasseB^2))^0.5)+centro(2);
end
yROI=real(yROI);

if options.display > 2
    hf.img_roi=figure();
    imagesc(immagine.img,immagine.bound),colormap(gray)
    hold on
    plot(xROI,yROI,'r')
    plot(centro(1),centro(2),'r+')
    title('AIF extraction - searching area')
end

% 2) DECIMAZIONE DEI VOXEL CANDIDATI
if options.display > 0
    disp('   Candidate voxel analysis')
end
% 2.1) Selezione a causa dell'area sotto la curva.
totCandidati=sum(sum(ROI));
totCandidatiDaTenere=ceil(totCandidati.*(1-pArea));
AUC=sum(AIFslice,3); % calcolo l'AUC di ogni voxel.
AUC=AUC.*ROI;
AUC(isinf(AUC))=0;

ciclo=true;
nCiclo=0;
AUCdown=min(min(AUC));
AUCup=max(max(AUC));
while ciclo
    nCiclo=nCiclo+1;
    soglia=0.5*(AUCup+AUCdown);
    nCandidati=sum(sum(AUC>soglia));
    
    if nCandidati==totCandidatiDaTenere
        ciclo=false;
    elseif nCandidati>totCandidatiDaTenere
        AUCdown=soglia;
    else
        AUCup=soglia;
    end
    if ((AUCup-AUCdown)<0.01)||(nCiclo>100)
        ciclo=false;
    end
end


ROIauc=2.*ROI-ROI.*(AUC>soglia); % Vale 2 per i voxel scartati, 1 per quelli tenuti.
if options.display > 2 
    disp(' ')
    disp(' Candidate voxel selection via AUC criteria')
    disp(['  Voxel initial amount: ' num2str(totCandidati)])
    disp(['  Survived voxels:      ' num2str(sum(sum(ROI)))])
    
    hf.aif=figure();
    subplot(2,3,1)
    imagesc(ROIauc)
    title('AUC criteria - survived voxel')
    set(gca,'xtick',[],'ytick',[])
    axis square
    
    subplot(2,3,4)
    plot(options.time,zeros(1,nT),'b-')
    hold on
    plot(options.time,zeros(1,nT),'r-')
    plot(options.time,zeros(1,nT),'k-')
    
    for c=1:nC
        for r=1:nR
            
            if ROIauc(r,c)==2
                plot(options.time,reshape(AIFslice(r,c,:),1,nT),'b-')
            end
        end
    end
    for c=1:nC
        for r=1:nR
            
            if ROIauc(r,c)==1
                plot(options.time,reshape(AIFslice(r,c,:),1,nT),'r-')
            end
        end
    end
    
    set(gca,...
        'FontSize',12)
    legend('Accepted','Reiected')
    title('AUC')
    xlabel('time')
end

ROI=ROI.*(AUC>soglia);



% 2.2) Selezione a causa del TTP
totCandidati=sum(sum(ROI));
totCandidatiDaTenere=ceil(totCandidati.*(1-pTTP));
[MC,TTP]=max(AIFslice,[],3);
TTP=TTP.*ROI;

ciclo=true;
soglia=1;
while ciclo
    if (sum(sum(TTP<soglia))-sum(sum(TTP==0)))>=totCandidatiDaTenere
        ciclo=false;
    else
        soglia=soglia+1;
    end
end


ROIttp=2*ROI-ROI.*(TTP<soglia); % Vale 2 per i voxel scartati, 1 per quelli tenuti.

if options.display > 2
    disp(' ')
    disp(' Candidate voxel selection via TTP criteria')
    disp(['  Voxel initial amount: ' num2str(totCandidati)])
    disp(['  Survived voxels:      ' num2str(sum(sum(ROI)))])
    
    figure(hf.aif);
    subplot(232)
    imagesc(ROIttp)
    title('TTP criteria - survived voxels')
    set(gca,'xtick',[],'ytick',[])
    axis square
    
    subplot(235)
    plot(options.time,zeros(1,nT),'b-')
    hold on
    plot(options.time,zeros(1,nT),'r-')
    plot(options.time,zeros(1,nT),'k-')
    for c=1:nC
        for r=1:nR
            
            if ROIttp(r,c)==1
                plot(options.time,reshape(AIFslice(r,c,:),1,nT),'r-')
            elseif ROIttp(r,c)==2
                plot(options.time,reshape(AIFslice(r,c,:),1,nT),'b-')
            end
        end
    end
    set(gca,...
        'FontSize',12)
    legend('Accepted','Reiected')
    title('TTP')
    xlabel('time')
end

ROI=ROI.*(TTP<soglia);


% 2.3) Seleziona in base all'indice di irregolarit�
totCandidati=sum(sum(ROI));
totCandidatiDaTenere=ceil(totCandidati.*(1-pReg));
REG=calcolaReg(AIFslice,options.time,ROI);

ciclo=true;
nCiclo=0;
REGdown=min(min(REG));
REGup=max(max(REG));
while ciclo
    nCiclo=nCiclo+1;
    soglia=0.5*(REGup+REGdown);
    nCandidati=sum(sum(REG>soglia));
    
    if nCandidati==totCandidatiDaTenere
        ciclo=false;
    elseif nCandidati<totCandidatiDaTenere
        REGup=soglia;
    else
        REGdown=soglia;
    end
    if ((REGup-REGdown)<0.001)||(nCiclo>=100)
        ciclo=false;
    end
end

ROIreg=2*ROI-ROI.*(REG>soglia); % Vale 2 per i voxel scartati, 1 per quelli tenuti.
if options.display > 2
    disp(' ')
    disp(' Candidate voxel selection via regularity criteria')
    disp(['  Voxel initial amount: ' num2str(totCandidati)])
    disp(['  Survived voxels: ' num2str(sum(sum(ROI)))])
    
    figure(hf.aif);
    subplot(2,3,3)
    imagesc(ROIreg)
    title('Regularity criteria - survived voxels')
    set(gca,'xtick',[],'ytick',[])
    axis square
    
    subplot(236)
    plot(options.time,zeros(1,nT),'b-')
    hold on
    plot(options.time,zeros(1,nT),'r-')
    plot(options.time,zeros(1,nT),'k-')
    for c=1:nC
        for r=1:nR
            
            if ROIreg(r,c)==1
                plot(options.time,reshape(AIFslice(r,c,:),1,nT),'r-')
            elseif ROIreg(r,c)==2
                plot(options.time,reshape(AIFslice(r,c,:),1,nT),'b-')
            end
        end
    end
    set(gca,...
        'FontSize',12)
    legend('Accepted','Reiected')
    title('I_r_e_g')
    xlabel('time')
end
ROI=ROI.*(REG>soglia);


if options.display > 2 
    hf.sel_voxel=figure();
    subplot(121)
    [posCok,posRok]=find(ROI);
    imagesc(immagine.img,immagine.bound),colormap(gray)
    hold on
    plot(xROI,yROI,'r')
    plot(centro(1),centro(2),'r+')
    plot(posRok,posCok,'r.','MarkerSize',1)
    title('Cadidate voxels')
    set(gca,'xtick',[],'ytick',[])
    axis square
    
    subplot(122)
    [posCno,posRno]=find(ROIiniziale-ROI);
    plot(options.time,reshape(AIFslice(posRok(1),posCok(1),:),1,nT),'r-')
    hold on
    plot(options.time,reshape(AIFslice(posRno(1),posCno(1),:),1,nT),'b-')
    for c=max([length(posCno),length(posCok)]):-1:2
        if c<=length(posCno)
            plot(options.time,reshape(AIFslice(posRno(c),posCno(c),:),1,nT),'b-')
        end
        if c<=length(posCok)
            plot(options.time,reshape(AIFslice(posRok(c),posCok(c),:),1,nT),'r-')
        end
    end
%     for c=1:length(posCno)
%         plot(options.time,reshape(AIFslice(posRno(c),posCno(c),:),1,nT),'b-')
%     end
%     for c=1:length(posCok)
%         plot(options.time,reshape(AIFslice(posRok(c),posCok(c),:),1,nT),'r-')
%     end
    legend('Accepted','Rejected')
    title('ROI voxels')
    xlabel('time')
end

% 3) APPLICAZIONE DELL'ALGORITMO DI CLUSTER PER LA RICERCA DELL'ARTERIALE
if options.display
    disp('   Arterial voxels extraction')
end

% 3.1) Preparazione della matrice contenente i dati
dati2D=zeros(sum(sum(ROI)),nT);
ind=find(ROI);
k=nR*nC;
for t=1:nT
    dati2D(:,t)=AIFslice(ind+k*(t-1));
end
maskAIF=ROI;


% 3.2) Applico l'algoritmo di Cluster Gerarchico in modo ricorsivo

ciclo=true;
nCiclo=0;
clear AIFslice

while ciclo
    nCiclo=nCiclo+1;
    if options.display > 2
        disp(' ')
        disp(' ------------------------------------')
        disp(['  CICLE N# ' num2str(nCiclo)])
    end
    
    % Applico il cluster gerarchico
    [vettCluster,centroidi]=clusterGerarchico(dati2D,2);
    
    % Confronto i cluster e scelgo quale tenere
    [MC1,TTP1]=max(centroidi(1,:));
    [MC2,TTP2]=max(centroidi(2,:));
    
    if (((max([MC1 MC2])-min([MC1 MC2]))/max([MC1 MC2]))<diffPicco)&&(TTP1~=TTP2)
        % La differenza tra i picchi � minore della soglia, scelgo
        % basandomi sul TTP
        clusterScelto=1+(TTP2<TTP1); % Il risultato vale 1 se TTP1>TTP2 e 2 se TTP2>TTP1
        
        if options.display > 2
            disp('  Cluster selected via TTP criteria')
            disp(['   Selected cluster: ' num2str(clusterScelto)])
        end
        
    else
        % Scelgo basandomi sulla differenza tra picchi
        clusterScelto=1+(MC2>MC1); % Il risultato vale 1 se MC1>MC2 e 2 se MC2>MC1
        
        if options.display > 2
            disp('  Cluster selected via MC criteria')
            disp(['   Selected cluster: ' num2str(clusterScelto)])
        end
    end
    
    if (sum(vettCluster==clusterScelto)<nVoxelMin)&&(sum(vettCluster==(3-clusterScelto))>=nVoxelMin)
        % La popolazione del cluster scelto � inferiore al numero minimo di
        % voxel accettati, mentre l'altro cluster ne ha a sufficenza.
        % Scelgo l'altro cluster.
        clusterScelto=3-clusterScelto; % Inverto il cluster scelto (se era 2 diventa 1, se era 1 diventa2)
        
        if options.display > 2
            disp('  Cluster selected switched because of minimum voxel bound')
            disp(['   Selected cluster: ' num2str(clusterScelto)])
        end
    end
    
    % Tengo solo i dati relativi al cluster scelto
    voxelScelti=(vettCluster==clusterScelto);
    indMask=find(maskAIF);
    maskAIF(indMask)=voxelScelti;
    
    indVoxel=find(voxelScelti);
    nL=length(indVoxel);
    dati2Dold=dati2D;
    dati2D=zeros(nL,nT);
    for t=1:nT
        dati2D(:,t)=dati2Dold(indVoxel,t);
    end
    
    if options.display > 2
        disp(' ')
        disp([' Resume cicle n# ' num2str(nCiclo)])
        disp(['  Voxel initial amount: ' num2str(length(indMask))])
        disp(['  Survived voxels:  ' num2str(nL)])
        disp(['  Cluster 1: MC       ' num2str(MC1)])
        disp(['             TTP      ' num2str(TTP1)])
        disp(['             voxel    ' num2str(sum(vettCluster==1))])
        disp(['  Cluster 2: MC       ' num2str(MC2)])
        disp(['             TTP      ' num2str(TTP2)])
        disp(['             voxel    ' num2str(sum(vettCluster==2))])
        disp(['  Selected cluster: ' num2str(clusterScelto)])
        
        
        eval(['hf.img_centr' num2str(nCiclo) '=figure();']);
        subplot(1,2,1)
        [posC,posR]=find(maskAIF);    
        imagesc(immagine.img,immagine.bound)
        colormap(gray)
        hold on
        plot(xROI,yROI,'r')
        plot(posR,posC,'r.','MarkerSize',1)
        title(['Cicle n#' num2str(nCiclo) ' - candidate voxels'])
        set(gca,'xtick',[],'ytick',[],'fontsize',12)
        axis square
        
        subplot(1,2,2)
        plot(options.time,centroidi,'k-')
        hold on
        plot(options.time,centroidi(clusterScelto,:),'r-')
        title('Cluster centroids')
        xlabel('time')
        set(gca,'fontsize',12)
        
    end
    
    
    
    % Controllo i criteri di uscita
    if (nL<=nVoxelMax)||(nCiclo>=100)
        ciclo=false;
    end
end

% 4) PREPARAZIONE DELL'UOTPUT

% 4.1) Salvo la ROI di ricerca
AIF.ROI.ind=find(ROI);
AIF.ROI.x=xROI;
AIF.ROI.y=yROI;

% 4.2) Salvo la posizione dei voxel scelti e la concentrazione media
AIFconc=centroidi(clusterScelto,:); % Campioni di concentrazione per l'AIF.
AIF.conc=AIFconc;
pos=1;
for r=1:nR
    for c=1:nC
        if maskAIF(r,c)==1
            AIF.voxels(pos,1)=r;
            AIF.voxels(pos,2)=c;
            pos=pos+1;
        end
    end
end


% 4.3) Calcolo il fit dell'arteriale con la gamma-variata (con ricircolo)
if options.display > 0
    disp('   Gamma variate fit computation')
end

pesi=0.01+exp(-AIFconc);              % Pesi per il calcolo del fit.

[MC TTP]=max(AIFconc);
pesi(TTP)=pesi(TTP)./10;
pesi(TTP-1)=pesi(TTP-1)./5;
pesi(TTP+1)=pesi(TTP+1)./2;


p  = {'t0' ''; 'alpha' '' ;'beta' ''; 'A' '' ; 'td' '';  'K' '' ; 'tao' '';'ExitFlag' ''};

[fitParameters_picco1,cv_est_parGV_picco1]=fitGV_picco1(AIFconc,pesi,options);
if options.aif.ricircolo
    [fitParameters_picco2,cv_est_parGV_picco2]=fitGV_picco2(AIFconc,pesi,fitParameters_picco1,options);
    fitParameters=[fitParameters_picco1(1:4) fitParameters_picco2(1:3)]';
    cv_est_parGV=[cv_est_parGV_picco1 cv_est_parGV_picco2];
else
    fitParameters=fitParameters_picco1(1:4);
    cv_est_parGV=cv_est_parGV_picco1;
end



AIF.fit.pesi=pesi;

AIF.fit.parameters=fitParameters;
AIF.fit.cv_est_parGV=cv_est_parGV;
if options.aif.ricircolo
    AIF.fit.gv=GVfunction(fitParameters,options);
else
    AIF.fit.gv=GVfunction_picco1(fitParameters,options);
end
AIF.fit.time=options.time;
AIF.conc=AIFconc;


%if options.display > 1

%    disp('Par: ' )

%    for j=1:size(p,1)
%        p{j,2}=AIF.fit.parameters(j);
%
%    end
%    disp(p);
%end

if options.display > 1    
    hf.fit_aif_final=figure();
    subplot(1,2,1);
    imagesc(real(immagine.img),immagine.bound),colormap(gray)
    hold on
    plot(AIF.ROI.x,AIF.ROI.y,'r','LineWidth',1)
    plot(AIF.voxels(:,2),AIF.voxels(:,1),'r.','MarkerSize',2)
    xlabel('','FontSize',10)
    legend('Searching area','Selected voxels')
    title('Arterial voxel location','FontSize',12)
    set(gca,'xtick',[],'ytick',[],'fontsize',10)  
    axis square
    
    subplot(1,2,2);
    
    plot(options.time,AIF.conc,'ko','MarkerSize',5)
    hold on
    plot(options.time,AIF.fit.gv,'k-','LineWidth',2)
    plot(options.time,dati2D,'r-')
    plot(options.time,AIF.conc,'ko','MarkerSize',5)
    plot(options.time,AIF.fit.gv,'k-','LineWidth',2)
    xlabel('[s]','FontSize',10)
    legend('AIF samples','GV function','Arterial voxels')
    title('AIF','FontSize',12)
    xlim([options.time(1) options.time(end)])
    hold off
    
end

if options.display > 2
    pause
    handles = fieldnames(hf);
    for i=1:size(handles,1)
        try
            eval(['close(hf.' handles{i,:} ');']);
        end
    end
end

%% ------------------------------------------------------------------------
function [REG]                         = calcolaReg(y,x,maschera)
% Calcola l'indice di irregolarit� dell'andamento della curva di
% concentrazione per ciascun voxel. L'indice viene calcolato normalizzando
% l'area a ! in modo da non penalizzare voxel con aree elevate.
% La formula usata per calcolare l'indice �:
%
% CTC=integrale((C"(t))^2 dt)



[nR,nC,nT]=size(y);
AUC=sum(y,3);
AUC=AUC+(AUC==1);
for t=1:nT
    y(:,:,t)=y(:,:,t)./AUC;
end

% calcolo della derivata seconda
derivata2=zeros(nR,nC,nT);
for r=1:nR
    for c=1:nC
        if maschera(r,c)==1
            for t=1:nT
                if (t>1)&&(t<nT)
                    % Caso standard
                    derivata2(r,c,t)=((y(r,c,t+1)-y(r,c,t))/(x(t+1)-x(t))-(y(r,c,t)-y(r,c,t-1))/(x(t)-x(t-1)))/(x(t)-x(t-1));
                elseif t==1
                    % Manca il campione precedente
                    derivata2(r,c,t)=(y(r,c,t+1)-y(r,c,t))/((x(t+1)-x(t))^2);
                else
                    % Manca il campione successivo
                    derivata2(r,c,t)=(y(r,c,t)-y(r,c,t-1))/((x(t)-x(t-1))^2);
                end
            end
        end
    end
end

% Calcolo dell'indice di irregolarit�
derivata2=derivata2.^2;
REG=trapz(x,derivata2,3);




%% ------------------------------------------------------------------------
function [vettAssegnazioni, centroidi] = clusterGerarchico(dati,nCluster)
% Applica l'algoritmo di cluster gerarchico ai dati e li divide in
% nCluster. Restutiusce un vettore contenente il numero del cluster cui �
% stato assegnato ciascun voxel e il centroide di tale cluster.

distance=pdist(dati);
tree=linkage(distance,'ward');
vettAssegnazioni=cluster(tree,'maxclust',nCluster);

nT=size(dati,2);
centroidi=zeros(nCluster,nT);
for k=1:nCluster
    ind=find(vettAssegnazioni==k);
    
    datiCluster=zeros(length(ind),nT);
    for t=1:nT
        datiCluster(:,t)=dati(ind,t);
    end
    
    centroidi(k,:)=mean(datiCluster,1);
end


%% ------------------------------------------------------------------------
function [GVparametri, cv_est_parGV]   = fitGV_picco1(dati,pesi,options_DSC)
% Calcola il fit del primo picco con una funzione gamma-variata.
% La funzione usata � descritta dalla formula:
%
% FP(t)=A*((t-t0)^alpha)*exp(-(t-t0)/beta)
%
% c(t)=FP(t)
%
% parametri: p=[t0 alpha beta A]
%
% L'ultimo parametro restituito rappresenta l'exitflag, che pu� assumere i
% seguenti valori:
%      1  LSQNONLIN converged to a solution X.
%      2  Change in X smaller than the specified tolerance.
%      3  Change in the residual smaller than the specified tolerance.
%      4  Magnitude search direction smaller than the specified tolerance.
%      5  Voxel nullo
%      0  Maximum number of function evaluations or of iterations reached.
%     -1  Algorithm terminated by the output function.
%     -2  Bounds are inconsistent.
%     -4  Line search cannot sufficiently decrease the residual along the
%         current search direction.

% OPZIONI STIMATORE
options             = optimset('lsqnonlin') ;
options.Display     = 'none'                ;
options.MaxFunEvals = 1000                 ;
options.MaxIter     = 1000                 ;
options.TolX        = 1e-4 ;
options.TolFun      = 1e-4 ;
%options.TolCon      = 1e-2 ;
%options.TolPCG      = 1e-8 ;
options.LargeScale  = 'on' ;
%options.DiffMinChange = 1e-18;
%options.DiffMaxChange = 2  ;

% STIME INIZIALI DEI PARAMETRI (modifica di DENIS)
% Alpha viene impostato a 5
alpha_init=5;

% t0 viene stimato sui dati iniziali. E' calcolato come l'ultimo istante in
% cui i dati rimangono in modulo inferiori al 5% del picco.
[MCdati,TTPpos]=max(dati);
TTPdati=options_DSC.time(TTPpos);
t0_init=options_DSC.time(find(dati(1:TTPpos)<=0.05*MCdati, 1, 'last' ));

% beta viene stimato sfruttando la relazione che TTP=t0+alpha*beta
beta_init=(TTPdati-t0_init)./alpha_init;

% Inizializzo i parametri [t0 alpha beta] e scelgo A in modo che la stima
% iniziale e i dati abbiano lo stesso massimo.
A_init= MCdati./max(GVfunction_picco1([t0_init; alpha_init; beta_init; 1],options_DSC));

% Valori iniziali dei parametri per la stima
% p  = [t0  alpha  beta  A]
p0   = [t0_init;   alpha_init;    beta_init;   A_init] ; % Valori iniziali
lb   = p0.*0.1; % Estremi inferiori
ub   = p0.*10 ; % Estremi superiori

if options_DSC.display>2
    h=figure();
    plot(options_DSC.time,dati,'ko',options_DSC.time,GVfunction_picco1(p0,options_DSC),'g-')
    title('First peak fit - initial values')
end

% Controllo i dati, devono essere vettori colonna
if size(options_DSC.time,1)==1
    % controlla che il vettore dei options.time sia un vettore colonna
    options_DSC.time=options_DSC.time';
end
if size(dati,1)==1
    % controlla che il vettore dei options.time sia un vettore colonna
    dati=dati';
end
if size(pesi,1)==1
    % controlla che il vettore dei options.time sia un vettore colonna
    pesi=pesi';
end

% MARCO
% AUMENTO LA PRECISIONE DEL PICCO
[MC TTP]=max(dati);
pesi(TTP)=pesi(TTP)./10;
pesi(TTP-1)=pesi(TTP-1)./2;


%TROVO FINE PRIMO PICCO (20% valore massimo)
i=TTP;
while dati(i)>0.2*dati(TTP)
    i=i+1;
end

%ADATTO I DATI PER "SOLO PRIMO PICCO"
dati_picco1=zeros(size(dati));
dati_picco1(1:i)=dati(1:i);

pesi_picco1=0.01+zeros(size(pesi));
pesi_picco1(1:i)=pesi(1:i);

% STIMATORE
ciclo=true;
nCiclo=0;
p=p0;
while ciclo
    nCiclo=nCiclo+1;
    [p, resNorm, residui, exitFlag,OUTPUT,LAMBDA,JACOBIAN] = lsqnonlin(@objFitGV_picco1, p, lb, ub, options, dati_picco1, pesi_picco1,options_DSC) ;
    
    if (nCiclo>=4)||(exitFlag>0)
        ciclo=false;
    end
end
GVparametri=p';

J=JACOBIAN;
covp=inv(J'*J);
var=diag(covp);
sd=sqrt(var);

cv_est_parGV=(sd./p*100)';

if options_DSC.display>2
    figure(h);
    hold on
    plot(options_DSC.time,GVfunction_picco1(p,options_DSC),'r-')
    title('First peak final fit')
    pause
    try
        close(h);
    end
end

%% ------------------------------------------------------------------------
function [out]                         = objFitGV_picco1(p,dati,pesi,options)
% Funzione obiettivo da minimizzare per la funzione fitGV_picco1
vett=GVfunction_picco1(p,options);

out=(vett-dati)./pesi;


%% ------------------------------------------------------------------------
function [GV]                          = GVfunction_picco1(p,options)
% Calcola la funzione gamma-variata definita dai parametri contenuti in p.
% La funzione gamma-variata � definita dalla formula:
%
% GV(t)=A*((t-t0)^alpha)*exp(-(t-t0)/beta)
%
% parametri: p=[t0 alpha beta A]

t0    = p(1);    % t0
alpha = p(2);    % alpha
beta  = p(3);    % beta
A     = p(4);    % A

nT=length(options.time);
GV=zeros(nT,1);
for cont=1:nT
    t=options.time(cont);
    if t>t0
        GV(cont)=A*((t-t0)^alpha)*exp(-(t-t0)/beta);
    end
end


%% ------------------------------------------------------------------------
function [GVparametri, cv_est_parGV]   = fitGV_picco2(dati,pesi,cost_picco1,options_DSC)
% Calcola il fit con una funzione gamma-variata.
% La funzione usata � descritta dalla formula:
%
% FP(t)=A*((t-t0)^alpha)*exp(-(t-t0)/beta)
%
% c(t)=FP(t) + FP(t-td) conv K*exp(-t/tao)
%
% parametri: p=[t0 alpha beta A td K tao]
%
% L'ultimo parametro restituito rappresenta l'exitflag, che pu� assumere i
% seguenti valori:
%      1  LSQNONLIN converged to a solution X.
%      2  Change in X smaller than the specified tolerance.
%      3  Change in the residual smaller than the specified tolerance.
%      4  Magnitude search direction smaller than the specified tolerance.
%      5  Voxel nullo
%      0  Maximum number of function evaluations or of iterations reached.
%     -1  Algorithm terminated by the output function.
%     -2  Bounds are inconsistent.
%     -4  Line search cannot sufficiently decrease the residual along the
%         current search direction.

% OPZIONI STIMATORE
options             = optimset('lsqnonlin') ;
options.Display     = 'none'                ;
options.MaxFunEvals = 10000                 ;
options.MaxIter     = 10000                 ;
options.TolX        = 1e-8 ;
options.TolFun      = 1e-8 ;
%options.TolCon      = 1e-2 ;
%options.TolPCG      = 1e-8 ;
options.LargeScale  = 'on' ;
%options.DiffMinChange = 1e-18;
%options.DiffMaxChange = 2  ;

% CONTROLLO DEI DATI
% Devono essere tutti vettori colonna
if size(options_DSC.time,1)==1
    % controlla che il vettore dei options.time sia un vettore colonna
    options_DSC.time=options_DSC.time';
end
if size(dati,1)==1
    % controlla che il vettore dei options.time sia un vettore colonna
    dati=dati';
end
if size(pesi,1)==1
    % controlla che il vettore dei options.time sia un vettore colonna
    pesi=pesi';
end
if size(cost_picco1,1)==1
    % controlla che il vettore dei options.time sia un vettore colonna
    cost_picco1=cost_picco1';
end

% PREPARO I DATI PER IL FIT
picco1=GVfunction_picco1(cost_picco1,options_DSC);
dati_picco2=dati-picco1; % I dati da fittare sono i residui del primo fit

%    pesi_picco2=0.01+exp(-dati_picco2); % Pesi per il calcolo del fit.

pesi_picco2=ones(length(dati_picco2),1); % Pesi per il calcolo del fit. (PESI UNIFORMI)
posTaglioPesi=min([find(dati>0.4*max(dati),1,'last'), 3+find(dati==max(dati))]);
% Riduco il peso dei dati prima del picco principale. Arrivo fino a quando
% la concentrazione non � scesa sotto il 40% del picco per evitare casi in
% cui il picco principale non fitta bene i dati successivi al picco e i
% residui potrebbero mostrare un picco fasullo.
pesi_picco2(1:posTaglioPesi)=1; 


% INIZIALIZZAZIONE PARAMETRI
% ricerca dei punti iniziali basata sui dati. Considero solo i dati
% dall'istante in cui le concentrazioni scendono sotto il 40% del picco per
% evitare residui troppo rumorosi non relativi al ricircolo. NB: il fit
% viene fatto con tutti i residui.
dati_x_stime_init=dati_picco2;
dati_x_stime_init(1:posTaglioPesi)=0;
dati_x_stime_init(find(dati_x_stime_init<0))=0;

% td_init viene calcolata come distanza tra l'istante del picco principale
% e la distanza del picco del ricircolo. Il picco del ricircolo viene
% individuato come picco dei dati meno la predizione del picco principale.
[maxPicco2,TTPpicco2]=max(dati_x_stime_init);
t0picco2 = find(dati_x_stime_init(1:TTPpicco2)<(0.1*max(dati_x_stime_init)),1,'last');
td_init = options_DSC.time(t0picco2)-cost_picco1(1);

% La stima iniziale di tao � fissata. 100 riesce a dare un ricircolo ben
% spalmato e che con i bound pu� diventare sia una dispersione nulla che
% portare ad una dispersione quasi completa.
tao_init=40;

% La stima di K viene fatta in modo che i massimi del ricircolo predetto e
% dei dati da fittare sia uguale.
ricircolo=GVfunction_ricircolo([cost_picco1; td_init; 1; tao_init],options_DSC);
K_init=max(dati_x_stime_init)./max(ricircolo);


% p  = [td   K   tao]
p = [td_init; K_init; tao_init ] ; % Valori iniziali


if options_DSC.display>2
    h=figure();
    plot(options_DSC.time,dati,'ko',options_DSC.time,picco1,'k-',options_DSC.time,GVfunction_ricircolo([cost_picco1; p],options_DSC),'g-')
    title('Recirculation fit - initial values')
end

% STIMATORE
ciclo=true;
nCiclo=0;
while ciclo
    ub=p.*[10; 10; 10];
    lb=p./[10; 10; 10];
    nCiclo=nCiclo+1;
    [p, resNorm, residui, exitFlag,OUTPUT,LAMBDA,JACOBIAN] = lsqnonlin(@objFitGV_picco2, p, lb, ub, options, dati_picco2,  pesi_picco2,cost_picco1,options_DSC) ;
    
    if (nCiclo>=4)||(exitFlag>0)
        ciclo=false;
    end
end
GVparametri=p';

J=JACOBIAN;
covp=inv(J'*J);
var=diag(covp);
sd=sqrt(var);

cv_est_parGV=(sd./p*100)';

if options_DSC.display>2
    figure(h);
    hold on
    plot(options_DSC.time,GVfunction_ricircolo([cost_picco1; p],options_DSC),'r-')
    title('Recirculation final fit')
    pause
    try
        close(h);
    end
end

%% ------------------------------------------------------------------------
function [out]                         = objFitGV_picco2(p,dati,pesi,cost_picco1,options)
% Funzione obiettivo da minimizzare per la funzione fitGV
vett=GVfunction_ricircolo([cost_picco1; p],options);

out=(vett-dati)./pesi;


%% ------------------------------------------------------------------------
function [ricircolo]                   = GVfunction_ricircolo(p,options)

% Calcola la funzione gamma-variata che descrive il ricircolo delle
% concentrazioni e definita dai parametri contenuti in p.
% La funzione gamma-variata � definita dalla formula:
%
% FP(t)= A*((t-t0)^alpha)*exp(-(t-t0)/beta)
%
% ricircolo(t)= FP(t-td) conv K*exp(-t/tao)
%
% I parametri sono passati nel sequente ordine:
% p= [t0 alpha beta A td K tao]
%
% NB: dato che la formula prevede una convoluzione, la griglia temporale
% lungo la quale viene calcolata la gamma-variata � molto pi� fitta della
% griglia finale.

t0    = p(1);    % t0
alpha = p(2);    % alpha
beta  = p(3);    % beta
A     = p(4);    % A
td    = p(5);    % td
K     = p(6);    % K
tao   = p(7);    % tao

% 1) Definizione della griglia virtuale necessaria per la convoluzione
TR=options.time(2)-options.time(1);
Tmax=max(options.time);
Tmin=min(options.time);
nT=length(options.time);

TRfine= TR/10;
tGrid=Tmin: TRfine : 2*Tmax;
nTfine=length(tGrid);

% 2) Calcolo le funzioni necessarie per il calcolo del ricircolo
% FP(t)        = A*((t-t0)^alpha)*exp(-(t-t0)/beta)
% disp(t)      = exp(-t/tao)
% ricircolo(t) = K * [FP(t-td) convoluto disp(t)]

% Inizializzazione dei vettori
picco2 = zeros(nTfine,1); % Picco del ricircolo
disp   = zeros(nTfine,1); % Dispersione del ricircolo

for cont=1:nTfine
    t=tGrid(cont);
    
    if t>t0+td
        % Calcolo di FP(t-td)
        picco2(cont)=K*((t-t0-td)^alpha)*exp(-(t-t0-td)/beta);
    end
    
    % Calcolo di disp(t)
    disp(cont)=exp(-t/tao);
end

% 3) Assemblo le componenti per ottenere la GV calcolata sulla griglia fine
ricircolo_fine=TRfine.*filter(picco2,1,disp);

% 4) Vado a campionare GV sugli istanti temporali richiesti in options.time
ricircolo=interp1(tGrid,ricircolo_fine,options.time);


%% ------------------------------------------------------------------------
function [GV]                          = GVfunction_picco2(p,options)

% Calcola le funzioni gamma-variata che descrivono l'andamento delle
% concentrazioni e definita dai parametri contenuti in p.
% La funzione complessiva � definita dalla formula:
%
% FP(t)= A*((t-t0)^alpha)*exp(-(t-t0)/beta)
%
% ricircolo(t)= FP(t-td) conv K*exp(-t/tao)
%
% I parametri sono passati nel sequente ordine:
% p= [t0 alpha beta A td K tao]
%
% NB: dato che la formula prevede una convoluzione, la griglia temporale
% lungo la quale viene calcolata la gamma-variata � molto pi� fitta della
% griglia finale.

FP=GVfunction_picco1(p(1:4),options);
ricircolo=GVfunction_ricircolo(p,options);

GV=FP+ricircolo;


%% ------------------------------------------------------------------------
function [GV]                          = GVfunction(p,options)
% Calcola la funzione gamma-variata definita dai parametri contenuti in p.
% La funzione gamma-variata � definita dalla formula:
%
% FP(t)=A*((t-t0)^alpha)*exp(-(t-t0)/beta)
%
% c(t)=FP(t) + FP(t-td) conv K*exp(-t/tao)
%
% parametri: p=[t0 alpha beta A td K tao]
% NB: dato che la formula prevede una convoluzione, la griglia temporale
% lungo la quale viene calcolata la gamma-variata � molto pi� fitta della
% griglia finale.

t0    = p(1);    % t0
alpha = p(2);    % alpha
beta  = p(3);    % beta
A     = p(4);    % A
td    = p(5);    % td
K     = p(6);    % K
tao   = p(7);    % tao

% 1) Definizione della griglia virtuale
TR=options.time(2)-options.time(1);
Tmax=max(options.time);
nT=length(options.time);

TRfine= TR/10;
tGrid=0: TRfine : 2*Tmax;
nTfine=length(tGrid);

% 2) Calcolo delle componenti di GV
% Divido la GV nelle sue componenti principali
picco1 = zeros(1,nTfine); % Picco principale
picco2 = zeros(1,nTfine); % Picco del ricircolo
disp   = zeros(1,nTfine); % Dispersione del ricircolo

for cont=1:nTfine
    t=tGrid(cont);
    
    if t>t0
        picco1(cont)=A*((t-t0)^alpha)*exp(-(t-t0)/beta);
    end
    
    if t>t0+td
        picco2(cont)=K*((t-t0-td)^alpha)*exp(-(t-t0-td)/beta);
    end
    
    disp(cont)=exp(-t/tao);
end

% 3) Assemblo le componenti per ottenere la GV calcolata sulla griglia fine
ricircolo=TRfine.*filter(picco2,1,disp);
conc=picco1+ricircolo;

% 4) Vado a campionare GV sui options.time richiesti
GV=zeros(1,nT);
for cont=1:nT
    [err,pos]=min(abs(tGrid-options.time(cont)));
    GV(cont)=conc(pos);
    
    if err>1
        disp('WARNING: approssimazione non buona.')
    end
end



%%