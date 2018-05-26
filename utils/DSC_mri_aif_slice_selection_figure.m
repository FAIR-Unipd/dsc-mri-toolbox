function varargout = DSC_mri_aif_slice_selection_figure(varargin)
% Visualizza l'immagine del volume cerebrale e chiede all'operatore di
% selezionare la slice sulla quale individuare la slice.

% DSC_MRI_AIF_SLICE_SELECTION_FIGURE M-file for DSC_mri_aif_slice_selection_figure.fig
%      DSC_MRI_AIF_SLICE_SELECTION_FIGURE, by itself, creates a new DSC_MRI_AIF_SLICE_SELECTION_FIGURE or raises the existing
%      singleton*.
%
%      H = DSC_MRI_AIF_SLICE_SELECTION_FIGURE returns the handle to a new DSC_MRI_AIF_SLICE_SELECTION_FIGURE or the handle to
%      the existing singleton*.
%
%      DSC_MRI_AIF_SLICE_SELECTION_FIGURE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DSC_MRI_AIF_SLICE_SELECTION_FIGURE.M with the given input arguments.
%
%      DSC_MRI_AIF_SLICE_SELECTION_FIGURE('Property','Value',...) creates a new DSC_MRI_AIF_SLICE_SELECTION_FIGURE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DSC_mri_aif_slice_selection_figure_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DSC_mri_aif_slice_selection_figure_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DSC_mri_aif_slice_selection_figure

% Last Modified by GUIDE v2.5 10-Jun-2010 11:18:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DSC_mri_aif_slice_selection_figure_OpeningFcn, ...
                   'gui_OutputFcn',  @DSC_mri_aif_slice_selection_figure_OutputFcn, ...
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



% --- Executes just before DSC_mri_aif_slice_selection_figure is made visible.
function DSC_mri_aif_slice_selection_figure_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DSC_mri_aif_slice_selection_figure (see VARARGIN)

volume=varargin{1};
options=varargin{2};

% Preparazione dei dati
handles.dati=sum(volume,4); % Immagine sommata per la visualizzazione
clear volume
[nR,nC,nS]=size(handles.dati);

% Individuo i bound ottimali per la visualizzazione delle slice
vettVolume=sort(handles.dati(find(handles.dati)));
handles.bound=[vettVolume(round(0.05*length(vettVolume))) vettVolume(round(0.9*length(vettVolume)))];
clear vettVolume

% Inizializzo le variabili di stato
handles.initialSlice=round(0.5*nS);

handles.nR=nR;
handles.nC=nC;
handles.nS=nS;

% Setup delle caratteristiche degli oggetti
 
% SLIDER
set(handles.slice_slider,'Max',handles.nS,'Min',1,'SliderStep',[1/(handles.nS-1) 1/(handles.nS-1)],'Value',handles.initialSlice);
% nSLICE_TEXT
set(handles.nSlice_text,'String',num2str(handles.initialSlice));
% SLICE_IMAGE

% image(handles.dati(:,:,handles.initialSlice),'Parent',handles.slice_image),colorbar
imagesc(real(handles.dati(:,:,handles.initialSlice)),'Parent',handles.slice_image)
caxis(handles.bound)
colormap('gray')
set(handles.slice_image,'Xtick',[],'Ytick',[])

% Choose default command line output for DSC_mri_aif_slice_selection_figure
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DSC_mri_aif_slice_selection_figure wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DSC_mri_aif_slice_selection_figure_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = get(handles.slice_slider,'Value');
close(handles.figure1)



% --- Executes on button press in ok_button.
function ok_button_Callback(hObject, eventdata, handles)
% hObject    handle to ok_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);


% --- Executes on slider movement.
function slice_slider_Callback(hObject, eventdata, handles)
% hObject    handle to slice_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

newSlice=round(get(hObject,'Value'));
set(hObject,'Value',newSlice);
set(handles.nSlice_text,'String',num2str(newSlice));
imagesc(real(handles.dati(:,:,newSlice)),'Parent',handles.slice_image)
caxis(handles.bound)
colormap('gray')
set(handles.slice_image,'Xtick',[],'Ytick',[])



% --- Executes during object creation, after setting all properties.
function slice_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slice_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function nSlice_text_Callback(hObject, eventdata, handles)
% hObject    handle to nSlice_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nSlice_text as text
%        str2double(get(hObject,'String')) returns contents of nSlice_text as a double


candidate=round(str2double(get(hObject,'String')));
% Se la stringa non contiene un numero valido allora str2double vale NaN
if (candidate>=1)&&(candidate<=handles.nS)
    newSlice=candidate;
else
    newSlice=get(handles.slice_slider,'Value');
end

set(handles.slice_slider,'Value',newSlice);
set(hObject,'String',num2str(newSlice));
imagesc(real(handles.dati(:,:,newSlice)),'Parent',handles.slice_image)
caxis(handles.bound)
colormap('gray')
set(handles.slice_image,'Xtick',[],'Ytick',[])

% --- Executes during object creation, after setting all properties.
function nSlice_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nSlice_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
