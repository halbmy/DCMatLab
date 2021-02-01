function varargout = bert2dpost(varargin)
% BERT2DPOST M-file for bert2dpost.fig
%      BERT2DPOST, by itself, creates a new BERT2DPOST or raises the existing
%      singleton*.
%
%      H = BERT2DPOST returns the handle to a new BERT2DPOST or the handle to
%      the existing singleton*.
%
%      BERT2DPOST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BERT2DPOST.M with the given input arguments.
%
%      BERT2DPOST('Property','Value',...) creates a new BERT2DPOST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before bert2dpost_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to bert2dpost_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help bert2dpost

% Last Modified by GUIDE v2.5 08-Jul-2008 17:35:08

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @bert2dpost_OpeningFcn, ...
                   'gui_OutputFcn',  @bert2dpost_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before bert2dpost is made visible.
function bert2dpost_OpeningFcn(hObject, eventdata, handles, varargin)
global filename
filename='*.zip;*.cfg';
if exist('post2d.last'), 
    fid=fopen('post2d.last','r');
    filename=fgetl(fid);
    fclose(fid);
end

% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to bert2dpost (see VARARGIN)

% Choose default command line output for bert2dpost
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes bert2dpost wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = bert2dpost_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in openfile.
function openfile_Callback(hObject, eventdata, handles)
global filename
[fname,pname]=uigetfile({'*.zip;*.cfg';'*.zip';'*.cfg'},'Choose result',filename);
if ~isstr(fname), return; end
filename=fullfile(pname,fname);
for i=1:4, if ishandle(i), close(i); end; end
bert2dpost('plotall_Callback',gcbo,[],guidata(gcbo))


% --- Executes on button press in plotall.
function plotall_Callback(hObject, eventdata, handles)
global Mesh N filename
[Mesh,N]=postmodel2d(filename);
if isfield(N,'response'),
    out=sprintf('%d iterations: RMS=%.1f%%, Chi^2=%.1f\n',Mesh.iter,rms(N.r,N.response),chi2(N.r,N.response,N.err,1));
else
    out=sprintf('%d iterations: fit unknown\n',Mesh.iter);
end
set(handles.status,'String',out);

% --- Executes on button press in cauto.
function cauto_Callback(hObject, eventdata, handles)
ison=get(hObject,'Value');
onoff='Off';
if ison==0, onoff='On'; end
set(handles.cmin,'Enable',onoff);
set(handles.cmax,'Enable',onoff);    


% --- Executes during object creation, after setting all properties.
function cmin_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function cmin_Callback(hObject, eventdata, handles)
%        str2double(get(hObject,'String')) returns contents of cmin as a double


% --- Executes during object creation, after setting all properties.
function cmax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function cmax_Callback(hObject, eventdata, handles)


% --- Executes on button press in ipcauto.
function ipcauto_Callback(hObject, eventdata, handles)
ison=get(hObject,'Value');
onoff='Off';
if ison==0, onoff='On'; end
set(handles.ipcmin,'Enable',onoff);
set(handles.ipcmax,'Enable',onoff);    


% --- Executes on button press in plotres.
function plotres_Callback(hObject, eventdata, handles)
global Mesh N
alfa=Mesh.alfa;
if get(handles.calfa,'Value')==0, alfa(:)=1; end
mal=struct('cauto',get(handles.cauto,'Value'),'cbar',get(handles.cbar,'Value'),...
    'cmin',str2double(get(handles.cmin,'String')),'cmax',str2double(get(handles.cmax,'String')),...
    'canot','Ohmm','oldstyle',1,'clog',get(handles.clog,'Value'));
if ~ishandle(1),
    figure(1);clf;set(1,'Renderer','zbuffer');
    set(gcf,'MenuBar','none','NumberTitle','off','Name','Inversion result');
    di=max(Mesh.node)-min(Mesh.node); 
    po=get(gcf,'Position');po(1)=20;po(2)=50;po(3)=800;po(4)=di(2)/di(1)*1.3*800+50;
    while po(4)>800, po(3:4)=po(3:4)/2; end
    set(gcf,'Position',po);
else
    figure(1);clf;
    set(gcf,'MenuBar','none','NumberTitle','off','Name','Inversion result');
end
tripatchmod(Mesh,Mesh.model,alfa,mal);
t=text(max(Mesh.node(:,1)),min(Mesh.node(:,2)),'BERT@resistivity.net');
set(t,'HorizontalAlignment','right','VerticalAlignment','bottom');


% --- Executes during object creation, after setting all properties.
function ipcmin_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ipcmin_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function ipcmax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function ipcmax_Callback(hObject, eventdata, handles)


% --- Executes on button press in plotphase.
function plotphase_Callback(hObject, eventdata, handles)
global Mesh N
alfa=Mesh.alfa;
if get(handles.calfa,'Value')==0, alfa(:)=1; end
mal=struct('cauto',get(handles.ipcauto,'Value'),'cbar',get(handles.cbar,'Value'),...
    'cmin',str2double(get(handles.ipcmin,'String')),'cmax',str2double(get(handles.ipcmax,'String')),...
    'canot','mrad','oldstyle',1);
if isfield(Mesh,'ipmodel'),
    figure(4);clf;
    set(gcf,'MenuBar','none','NumberTitle','off','Name','IP model');
    tripatchmod(Mesh,Mesh.ipmodel,alfa,mal);
    t=text(max(Mesh.node(:,1)),min(Mesh.node(:,2)),'BERT@resistivity.net');
    set(t,'HorizontalAlignment','right','VerticalAlignment','bottom');
end


% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)
global filename
[pname,basename,ext]=fileparts(filename);
if ishandle(1), epsprint(1,[pname filesep basename '-model'],1); end
if ishandle(2), epsprint(2,[pname filesep basename '-data'],1); end
if ishandle(3), epsprint(3,[pname filesep basename '-topoeff'],1); end
if ishandle(4), epsprint(4,[pname filesep basename '-ipmodel'],1); end


% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
global filename
for i=1:4, if ishandle(i), close(i); end; end
fid=fopen('post2d.last','w');
fprintf(fid,'%s\n',filename);
fclose(fid);
delete(gcbf);


% --- Executes on button press in calfa.
function calfa_Callback(hObject, eventdata, handles)


% --- Executes on button press in cbar.
function cbar_Callback(hObject, eventdata, handles)


% --- Executes on button press in clog.
function clog_Callback(hObject, eventdata, handles)
% hObject    handle to clog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of clog


