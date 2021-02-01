function varargout = m_mod(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @m_mod_OpeningFcn, ...
                   'gui_OutputFcn',  @m_mod_OutputFcn, ...
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

% --- Executes just before m_mod is made visible.
function m_mod_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to m_mod (see VARARGIN)

% Choose default command line output for m_mod
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
iconify(hObject);

if strcmp(get(hObject,'Visible'),'off')
    global Model N Mod
    Mod=Model;
    init_gui(hObject, handles);
end

% UIWAIT makes m_mod wait for user response (see UIRESUME)
%uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = m_mod_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function init_gui(fig_handle, handles)
global Mod
nlay=length(Mod.z)-1;
set(handles.nlay,'String',num2str(nlay));
if isfield(Mod,'x0'), set(handles.xnull,'String',num2str(Mod.x0));
else set(handles.xnull,'String',num2str(min(Mod.x))); end
if isfield(Mod,'Y0'), set(handles.ynull,'String',num2str(Mod.y0));
else set(handles.ynull,'String',num2str(min(Mod.y))); end
set(handles.dx,'String',num2str(Mod.dx));
set(handles.dy,'String',num2str(Mod.dy));
set(handles.z,'String',sprintf('%.1f ',Mod.z));
set(handles.bg,'String',sprintf('%d ',Mod.Bg));
set(handles.ncells,'String',[num2str(Mod.ncells) ' cells']),
set(handles.nx,'String',sprintf('%d ',Mod.nx));
set(handles.ny,'String',sprintf('%d ',Mod.ny));

% --- Executes during object creation, after setting all properties.
function nlay_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function nlay_Callback(hObject, eventdata, handles)
global Mod N
nlay=round(str2num(get(handles.nlay,'String')));
if (nlay>0)&&(nlay~=(length(Mod.z)-1)),
    dx=str2num(get(handles.dx,'String'));
    dy=str2num(get(handles.dy,'String'));
    Mod=create3dmod(N,nlay,dx,dy);
    init_gui(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function z_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function z_Callback(hObject, eventdata, handles)
global Mod
z=str2num(get(handles.z,'String'));
if (length(z)==length(Mod.z)),
    Mod.z=z;
end
set(handles.bg,'String',sprintf('%d ',round(Mod.Bg)));

% --- Executes during object creation, after setting all properties.
function dx_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function dx_Callback(hObject, eventdata, handles)
global Mod Model N
dx=str2num(get(handles.dx,'String'));
dy=str2num(get(handles.dy,'String'));
nlay=round(str2num(get(handles.nlay,'String')));
if isempty(dx),
    set(handles.dx,'String',str2num(Model.dx));
else
    %     Mod.dx=dx;
    Mod=create3dmod(N,nlay,dx,dy);
    init_gui(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function dy_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function dy_Callback(hObject, eventdata, handles)
global Mod Model N
dx=str2num(get(handles.dx,'String'));
dy=str2num(get(handles.dy,'String'));
nlay=round(str2num(get(handles.nlay,'String')));
if isempty(dy),
    set(handles.dy,'String',num2str(Model.dy));
else
    %     Mod.dy=dy;
    Mod=create3dmod(N,nlay,dx,dy);
    init_gui(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function bg_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function bg_Callback(hObject, eventdata, handles)
global Mod
Bg=str2num(get(handles.bg,'String'));
if (length(Bg)==length(Mod.Bg)),
    Mod.Bg=round(Bg);
    for k=1:length(Mod.Bg),
        Mod.M{k}(:)=Mod.Bg(k);
    end
end
set(handles.bg,'String',sprintf('%d ',Mod.Bg));

% --- Executes during object creation, after setting all properties.
function xnull_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function xnull_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ynull_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function ynull_Callback(hObject, eventdata, handles)

% --- Executes on button press in show.
function show_Callback(hObject, eventdata, handles)
global Mod
figure(1);clf;
draw3dmodel(Mod,struct('log',1,'nu',0,'nv',0));

% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
if ishandle(11), delete(11); end
delete(gcbf);

% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
global Model Mod
Model=Mod;
clear global Mod
if ishandle(11), delete(11); end
delete(gcbf);

% --- Executes during object creation, after setting all properties.
function nx_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function nx_Callback(hObject, eventdata, handles)
global Model Mod
nx=round(str2num(get(handles.nx,'String')));
if length(nx)==length(Mod.nx),
    Mod.nx=nx;
    for k=1:length(Mod.M),
      nnx=fix(size(Model.M{1},1)*Model.nx(1)/Mod.nx(k));
      Mod.M{k}=ones(nnx,size(Mod.M{k},2))*Mod.Bg(k);
  end
end

% --- Executes during object creation, after setting all properties.
function ny_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function ny_Callback(hObject, eventdata, handles)
global Model Mod
ny=round(str2num(get(handles.ny,'String')));
if length(ny)==length(Mod.ny),
    Mod.ny=ny;
    for k=1:length(Mod.M),
      nny=fix(size(Model.M{1},2)*Model.ny(1)/Mod.ny(k));
      Mod.M{k}=ones(size(Mod.M{k},1),nny)*Mod.Bg(k);
  end
end