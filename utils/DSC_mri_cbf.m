function [cbf]=DSC_mri_cbf(conc,aif,mask,options)
%Funzione del pacchetto DSC_mri - DSC_mri_cbf
%Autore: Castellaro Marco - Università di Padova - DEI
%
%Calcola le mappe parametriche di Cerebral Blood Flow (CBF) per un soggetto
%
%Parametri in ingresso - conc (Matrice 4D), contiene gli andamenti delle
%                        concentrazioni DSC di tutti i voxel. 
%                      - aif, andamento delle concentrazioni nel sito
%                        scelto come arteriale
%                      - mask (Matrice 3D), contiene la matrice utilizzata
%                        per mascherare il volume cerebrale da analizzare
%
%Options è la sruct che contiene le opzioni del metodo, quelli
%significativi sono:
%
%Il calcolo del parametro CBF richiede di effettuare una operazione di
%deconvoluzione, nel toolbox sono previste alcune metodologie adatte a
%questo scopo: SVD e cSVD (vers. Block-circulant).
%
%options.deconv.method - Deve essere un vettore di celle che contiene i
%                        nomi degli algoritmi con cui si intende realizzare
%                        l'analisi. Es. SVD  cSVD o SS
%
%                       Per ogni metodo deve essere inserita nelle opzioni
%                       una struct con i parametri caratteristici del
%                       metodo stesso, come in questo esempio:
%
%                       options.deconv.<nome metodo>.<par_1>
%                       options.deconv.<nome metodo>.<...>
%                       options.deconv.<nome metodo>.<par_n>
%
%options.time - Rappresenta il vettore dei tempi dell'esame DSC (ogni
%               campione rappresenta l'acquisizione di un intero volume cerebrale
%
%options.display - livello 1 Mostra l'avanzamento dell'elaborazione, 
%                  livello 2 Da anche informazioni sui parametri settati
%                  per gli algoritmi utilizzati
%
%Parametri in uscita: 
%cbf - diverse sub-struct, una per ogni metodo che si è scelto di utilizzare, 
%      ogni sub-struct contiene un campo map che contraddistingue la mappa 
%      di cbv calcolata, ad esempio:
%      cbf.<nome metodo>.map
%
%      residual, matrice 4D che contiene i residui (è un campo opzionale che
%      va richiesto nelle opzioni: parametro options.deconv.<nome metodo>.residual) 
%      ad esempio:
%      cbf.<nome metodo>.residual


% E' possibile aggiungere un metodo di deconvoluzione, per poterlo fare è 
% necessario aggiornare la variabile method, aggiungendo la stringa 
% identificativa e un caso per la chiamata del nuovo metodo, inoltre va
% scritta una funzione che realizzi il metodo e chiamata DSC_mri_<nome_metodo>

%%
%AGGIORNARE SE SI INTENDE AGGIUNGERE UN NUOVO METODO
method={'SVD';'cSVD';'oSVD'};

if options.display > 0
    disp('   CBF');
end


for alg=1:size(options.deconv.method,1)
    switch options.deconv.method{alg,:}
        case method{1,:} %SVD
            cbf.svd=DSC_mri_SVD(conc,aif,mask.data,options);
        case method{2,:} %cSVD
            cbf.csvd=DSC_mri_cSVD(conc,aif,mask.data,options);
        case method{3,:} %oSVD          
            cbf.osvd=DSC_mri_oSVD(conc,aif,mask.data,options);
            
        %AGGIORNARE SE SI INTENDE AGGIUNGERE UN NUOVO METODO
        %case method{n,:} 
        %cbf.<nome metodo>=DSC_mri_<nome_metodo>(...);
        
        otherwise %metodo non riconosciuto
            str_errore='';
            for m=1:size(method,1)
                %DA SISTEMARE
                str_errore=strcat(str_errore,method{m,:},' or ');
            end
            error(['Deconvolution method not recognized. Use: ' str_errore])
    end
end