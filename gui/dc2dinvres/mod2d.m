function varargout = mod2d(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mod2d_OpeningFcn, ...
                   'gui_OutputFcn',  @mod2d_OutputFcn, ...
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


% --- Executes just before mod2d is made visible.
function mod2d_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mod2d (see VARARGIN)

% Choose default command line output for mod2d
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
iconify(hObject);
global N Mod MAL FOR
load setup
Mod.x=0:10;
Mod.z=0:5;
Mod.Lay=100;
Mod.M=ones(length(Mod.x)-1,length(Mod.z)-1)*Mod.Lay(1);
Mod.M=[];Mod.Lay=0;
N.elec=[3 0; 7 0];
N.a=1;N.m=2;N.n=0;N.b=0;N.r=100;N.k=getkonf(N);

% UIWAIT makes mod2d wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mod2d_OutputFcn(hObject, eventdata, handles)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in forward.
function forward_Callback(hObject, eventdata, handles)
global N Mod FOR
N.r=dcfwd2d(Mod.x,Mod.z,Mod.M,Mod.Lay,N,FOR);
mod2d('showdata_Callback',gcbo,[],guidata(gcbo));

% --- Executes on button press in createdata.
function createdata_Callback(hObject, eventdata, handles)
global Mod N
uiwait(ncreate);
mod2d('showdata_Callback',gcbo,[],guidata(gcbo));
if isempty(Mod.M),
    [Mod.x,Mod.z,Mod.M]=modelfromdata2d(N);
    mod2d('showmodel_Callback',gcbo,[],guidata(gcbo));
end


function showdata_Callback(hObject, eventdata, handles)
global N
axes(handles.data(end));
showdata2d(N);
%im=get(handles.data,'Children')
%get(im(2))

function showmodel_Callback(hObject, eventdata, handles)
global N Mod MAL
axes(handles.model(end));
draw2dmodel(Mod.x,Mod.z,Mod.M,MAL);
im=get(handles.model(end),'Children');
for i=1:length(im),
    ty=get(im(i),'Type');
    if strcmp(ty,'image')|strcmp(ty,'surface'),
        set(im(i),'ButtonDownFcn','mod2d(''model_ButtonDownFcn'',gcbo,[],guidata(gcbo))')
    end
end

% --- Executes on mouse press over axes background.
function model_ButtonDownFcn(hObject, eventdata, handles)
global Mod
cp=get(handles.model(end),'CurrentPoint');
nx=max(find(Mod.x<cp(1,1)));
nz=max(find(Mod.z<cp(1,2)));
if (nx>=1)&(nx<length(Mod.x))&(nz>=1)&(nz<length(Mod.z)), %found
    nu=str2num(get(handles.widerstand,'String'));
    if isempty('nu'), nu=0; end
    if nu>0,
        Mod.M(nx,nz)=nu;
        mod2d('showmodel_Callback',gcbo,[],guidata(gcbo));
        set(handles.widerstand,'String',num2str(nu));
    else
        msgbox(sprintf('rho=%.1f (Mod.x=%.1f-%.1f, Mod.z=%.1f-%.1f)\n',Mod.M(nx,nz),Mod.x(nx),Mod.x(nx+1),Mod.z(nz),Mod.z(nz+1)),...
            'Resistivity cell');
    end    
end


% --- Executes during object creation, after setting all properties.
function widerstand_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function widerstand_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function options_Callback(hObject, eventdata, handles)
% --------------------------------------------------------------------
function optmal_Callback(hObject, eventdata, handles)
uiwait(s_mal);
mod2d('showmodel_Callback',gcbo,[],guidata(gcbo));

% --------------------------------------------------------------------
function optfor_Callback(hObject, eventdata, handles)
uiwait(s_for);

% --------------------------------------------------------------------
function exit_Callback(hObject, eventdata, handles)
delete(gcbf);


% --------------------------------------------------------------------
function noisifydata_Callback(hObject, eventdata, handles)
global N
daterr;
noise=randn(size(N.r)).*N.err;
N.r=N.r.*(1+noise);
mod2d('showdata_Callback',gcbo,[],guidata(gcbo));


% --------------------------------------------------------------------
function editmodel_Callback(hObject, eventdata, handles)
s_mal;
mod2d('showmodel_Callback',gcbo,[],guidata(gcbo));

% --------------------------------------------------------------------
function createmodel_Callback(hObject, eventdata, handles)
global Mod N
[Mod.x,Mod.z,Mod.M]=modelfromdata2d(N);
mod2d('showmodel_Callback',gcbo,[],guidata(gcbo));
