function varargout = DSC_mri_show_results(varargin)
%DSC_mri_show_results
%ultima modifica: Denis Peruzzo 11/06/2010
%Funzione del pacchetto DSC_mri
%Autore: Denis Peruzzo - Università di Padova - DEI
%
%Funzione che permette di visualizzare i risultati ottenuti con il
%pacchetto DSC_mri. Istruzione di chiamata standard:
%
% DSC_mri_show_results(CBVmap,CBFmap,MTTmap,TTPmap,mask,aif,conc,s0);
%
%Parametri in ingresso: 
% - CBVmap: matrice contenente le mappe di CBV.
% - CBFmap: variabile strutturata contenente le informazioni, le mappe e le
%           funzioni residuo ottenute dal calcolo del CBF.
% - MTTmap: variabile strutturata contenente tutte le mappe di MTT.
% - TTPmap: matrice contenente le mappe di TTP.
% - mask:   variabile strutturata contenente le informazioni relative alle
%           maschere utilizzate.
% - aif:    variabile strutturata contenente tutte le informazioni relative
%           alla funzione di ingresso arteriale.
% - conc:   matrice contenente l'andamento della concentrazione di
%           tracciante per ciascun voxel.
% - s0:     matrice contenente il valore di S0 di ciascun voxel.
%
% NB: tutte le variabili devono avere la stessa struttura fornita in uscita
% dal software DSC_mri

% DSC_MRI_SHOW_RESULTS M-file for DSC_mri_show_results.fig
%      DSC_MRI_SHOW_RESULTS, by itself, creates a new DSC_MRI_SHOW_RESULTS or raises the existing
%      singleton*.
%
%      H = DSC_MRI_SHOW_RESULTS returns the handle to a new DSC_MRI_SHOW_RESULTS or the handle to
%      the existing singleton*.
%
%      DSC_MRI_SHOW_RESULTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSC_MRI_SHOW_RESULTS.M with the given input arguments.
%
%      DSC_MRI_SHOW_RESULTS('Property','Value',...) creates a new DSC_MRI_SHOW_RESULTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DSC_mri_show_results_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DSC_mri_show_results_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DSC_mri_show_results

% Last Modified by GUIDE v2.5 11-Jun-2010 12:33:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DSC_mri_show_results_OpeningFcn, ...
                   'gui_OutputFcn',  @DSC_mri_show_results_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%% ------------------------------------------------------------------------
%                figure1: 241.0046
%      select_voxel_text: 71.0050
%            voxel_table: 70.0050
%           MTTmap_title: 69.0050
%           CBVmap_title: 68.0050
%           CBFmap_title: 67.0050
%        data_plot_title: 66.0050
%     residue_plot_title: 65.0050
%                  text2: 64.0050
%           voxel_button: 63.0050
%           figure_title: 62.0050
%     deconvolution_list: 61.0050
%            exit_button: 60.0050
%             slice_edit: 59.0051
%           slice_slider: 267.0046
%              data_plot: 262.0046
%           residue_plot: 257.0046
%                 MTTmap: 252.0046
%                 CBVmap: 247.0046
%                 CBFmap: 242.0046
%                 output: 241.0046
%%


% --- Executes just before DSC_mri_show_results is made visible.
function DSC_mri_show_results_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DSC_mri_show_results (see VARARGIN)

% PREPARAZIONE DEI DATI INIZIALI PASSATI IN INGRESSO
% ordine dei parametri in ingresso:
% [CBVmap,CBFmap,MTTmap,TTPmap,mask,aif,conc,s0]
handles.dati.CBV  = varargin{1};
handles.dati.CBF  = varargin{2};
handles.dati.MTT  = varargin{3};
handles.dati.TTP  = varargin{4};
handles.dati.mask = varargin{5}.data;
handles.dati.aif  = varargin{6};
handles.dati.conc = varargin{7};
[nR,nC,nS,nT]=size(handles.dati.conc);
handles.dati.nR=nR;
handles.dati.nC=nC;
handles.dati.nS=nS;
handles.dati.nT=nT;




% SETTING DEI VALORI INIZIALI DI SLIDER ED ELENCHI
handles.startingSlice=round(0.5*nS);
handles.startingVoxel=[round(0.5*nR) round(0.5*nC)];

% VALORI INIZIALI PARAMETRI DI VISUALIZZAZIONE
% slice slider
set(handles.slice_slider,'Max',nS,'Min',1,'SliderStep',[1/(nS-1) 1/(nS-1)],'Value',handles.startingSlice);

% slice text
set(handles.slice_edit,'String',num2str(handles.startingSlice));

% deconvolution list
decList=fields(handles.dati.CBF);
set(handles.deconvolution_list,'String',decList,'Value',1);

% select voxel text
set(handles.select_voxel_text,'String',[num2str(handles.startingVoxel(1)) '-' num2str(handles.startingVoxel(2))]);

% voxel table
voxTable{1,1}=handles.dati.CBF.svd.map(handles.startingVoxel(1),handles.startingVoxel(2),handles.startingSlice); %CBF
voxTable{2,1}=handles.dati.CBV(handles.startingVoxel(1),handles.startingVoxel(2),handles.startingSlice); %CBV
voxTable{3,1}=handles.dati.MTT.svd(handles.startingVoxel(1),handles.startingVoxel(2),handles.startingSlice); %MTT
voxTable{4,1}=handles.dati.TTP(handles.startingVoxel(1),handles.startingVoxel(2),handles.startingSlice); %TTP
set(handles.voxel_table,'Data',voxTable);

aggiornaPlot(handles);

% Choose default command line output for DSC_mri_show_results
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DSC_mri_show_results wait for user response (see UIRESUME)
uiwait(handles.figure1);
%% ------------------------------------------------------------------------

function []=aggiornaPlot(handles)
% Aggiorna tutti i plot della figura

% LEGGO LA POSIZIONE DI COSA VISUALIZZARE
slice=get(handles.slice_slider,'Value');
stringaVoxel=get(handles.select_voxel_text,'String');
row=str2double(stringaVoxel(1:find(stringaVoxel=='-')-1));
col=str2double(stringaVoxel(find(stringaVoxel=='-')+1:length(stringaVoxel)));
decList=get(handles.deconvolution_list,'String');
dec=decList{get(handles.deconvolution_list,'Value'),1};
vettInd= find(handles.dati.mask);

% MAPPE DI PERFUSIONE
% estraggo le mappe
eval(['mappaCBF=handles.dati.CBF.' dec '.map;']);
eval(['mappaCBV=handles.dati.CBV;']);
eval(['mappaMTT=handles.dati.MTT.' dec ';']);
eval(['mappaTTP=handles.dati.TTP;']);

% CBF
% Calcolo i bound
vettCBF=sort(mappaCBF(vettInd));
CBFbound=[0 vettCBF(round(0.95*length(vettCBF)))];

mappa=mappaCBF;
mappa(find(mappaCBF>CBFbound(2)))=CBFbound(2);
imagesc(mappa(:,:,slice),'Parent',handles.CBFmap)
colormap('hot')
hold(handles.CBFmap,'on')
plot(handles.CBFmap,col,row,'go','MarkerSize',5,'LineWidth',2)
set(handles.CBFmap,'Xtick',[],'Ytick',[])

% CBV
% Calcolo i bound
vettCBV=sort(mappaCBV(vettInd));
CBVbound=[0 vettCBV(round(0.95*length(vettCBV)))];

mappa=mappaCBV;
mappa(find(mappaCBV>CBVbound(2)))=CBVbound(2);
imagesc(mappa(:,:,slice),'Parent',handles.CBVmap)
colormap('hot')
hold(handles.CBVmap,'on')
plot(handles.CBVmap,col,row,'go','MarkerSize',5,'LineWidth',2)
set(handles.CBVmap,'Xtick',[],'Ytick',[])

% MAPPE DI PERFUSIONE
% MTT
% Calcolo i bound
vettMTT=sort(mappaMTT(vettInd));
MTTbound=[0 vettMTT(round(0.95*length(vettMTT)))];

mappa=mappaMTT;
mappa(find(mappaMTT>MTTbound(2)))=MTTbound(2);
imagesc(mappa(:,:,slice),'Parent',handles.MTTmap)%,'Clim',MTTbound)
colormap('hot')
hold(handles.MTTmap,'on')
plot(handles.MTTmap,col,row,'go','MarkerSize',5,'LineWidth',2)
set(handles.MTTmap,'Xtick',[],'Ytick',[])

% PLOT DEL VOXEL
tempi=handles.dati.aif.fit.time;
try
    eval(['residui=handles.dati.CBF.' dec '.residual;'])
    vettRes=reshape(residui(row,col,slice,1:handles.dati.nT),1,handles.dati.nT);
    
catch
    vettRes=zeros(1,handles.dati.nT);
end
vettAIF=handles.dati.aif.fit.gv;
plot(handles.residue_plot,tempi,vettRes,'r+-')
vettConc=reshape(handles.dati.conc(row,col,slice,1:handles.dati.nT),1,handles.dati.nT);
vettRiconv=filter(vettAIF,1,vettRes).*(tempi(2)-tempi(1));
plot(handles.data_plot,tempi,vettConc,'ko',tempi,vettRiconv,'r-')

% STATISTICHE RELATIVE AL VOXEL
voxTable{1,1}=mappaCBF(row,col,slice); %CBF
voxTable{2,1}=mappaCBV(row,col,slice); %CBV
voxTable{3,1}=mappaMTT(row,col,slice); %MTT
voxTable{4,1}=mappaTTP(row,col,slice); %TTP
set(handles.voxel_table,'Data',voxTable);
    





% --- Outputs from this function are returned to the command line.
function varargout = DSC_mri_show_results_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
close(handles.figure1)
%% ------------------------------------------------------------------------


% --- Executes on slider movement.
function slice_slider_Callback(hObject, eventdata, handles)
% hObject    handle to slice_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
newSlice=round(get(hObject,'Value'));
set(hObject,'Value',newSlice);
set(handles.slice_edit,'String',num2str(newSlice));

aggiornaPlot(handles);
%% ------------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function slice_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slice_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
%% ------------------------------------------------------------------------


function slice_edit_Callback(hObject, eventdata, handles)
% hObject    handle to slice_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slice_edit as text
%        str2double(get(hObject,'String')) returns contents of slice_edit as a double
candidate=round(str2double(get(hObject,'String')));
% Se la stringa non contiene un numero valido allora str2double vale NaN
if (candidate>=1)&&(candidate<=handles.dati.nS)
    newSlice=candidate;
else
    newSlice=get(handles.slice_slider,'Value');
end

set(handles.slice_slider,'Value',newSlice);
set(handles.slice_edit,'String',num2str(newSlice));

aggiornaPlot(handles);
%% ------------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function slice_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slice_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% ------------------------------------------------------------------------


% --- Executes on button press in exit_button.
function exit_button_Callback(hObject, eventdata, handles)
% hObject    handle to exit_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
%% ------------------------------------------------------------------------


% --- Executes on selection change in deconvolution_list.
function deconvolution_list_Callback(hObject, eventdata, handles)
% hObject    handle to deconvolution_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns deconvolution_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from deconvolution_list
aggiornaPlot(handles);
%% ------------------------------------------------------------------------


% --- Executes during object creation, after setting all properties.
function deconvolution_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deconvolution_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in voxel_button.
function voxel_button_Callback(hObject, eventdata, handles)
% hObject    handle to voxel_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[x,y]=ginput(1);
cand_row=round(y);
cand_col=round(x);
if (cand_col>0)&&(cand_col<=handles.dati.nC)&&(cand_row>0)&&(cand_row<=handles.dati.nR)
    set(handles.select_voxel_text,'String',[num2str(cand_row) '-' num2str(cand_col)]);
end
aggiornaPlot(handles);
%% ------------------------------------------------------------------------


% --- Executes when selected cell(s) is changed in voxel_table.
function voxel_table_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to voxel_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
%% ------------------------------------------------------------------------
