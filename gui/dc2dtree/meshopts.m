function varargout = meshopts(varargin)
% MESHOPTS M-file for meshopts.fig
%      MESHOPTS, by itself, creates a new MESHOPTS or raises the existing
%      singleton*.
%
%      H = MESHOPTS returns the handle to a new MESHOPTS or the handle to
%      the existing singleton*.
%
%      MESHOPTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MESHOPTS.M with the given input arguments.
%
%      MESHOPTS('Property','Value',...) creates a new MESHOPTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before meshopts_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to meshopts_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help meshopts

% Last Modified by GUIDE v2.5 03-Dec-2005 11:23:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @meshopts_OpeningFcn, ...
                   'gui_OutputFcn',  @meshopts_OutputFcn, ...
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


% --- Executes just before meshopts is made visible.
function meshopts_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to meshopts (see VARARGIN)

% Choose default command line output for meshopts
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global Meshopts
if isfield(Meshopts,'parasmooth'), 
    set(handles.pasmooth,'Value',Meshopts.parasmooth>0); 
else
    set(handles.pasmooth,'Value',0); 
end
if isfield(Meshopts,'primsmooth'), 
    set(handles.prsmooth,'Value',Meshopts.primsmooth>0); 
else
    set(handles.prsmooth,'Value',0); 
end
if isfield(Meshopts,'paraquality'),
    set(handles.paautoquality,'Value',1);
    set(handles.paquality,'Enable','on','String',num2str(Meshopts.paraquality));
else
    set(handles.paquality,'String','33.5');
end
if isfield(Meshopts,'primquality'),
    set(handles.prautoquality,'Value',1);
    set(handles.prquality,'Enable','on','String',num2str(Meshopts.primquality));
else
    set(handles.prquality,'String','33.5');
end
if isfield(Meshopts,'pararefine'),
    set(handles.paautorefine,'Value',1);
    set(handles.parefine,'Enable','on','String',num2str(Meshopts.pararefine));
else
    set(handles.parefine,'String','0.2');
end
if isfield(Meshopts,'primrefine'),
    set(handles.prautorefine,'Value',1);
    set(handles.prrefine,'Enable','on','String',num2str(Meshopts.primrefine));
else
    set(handles.prrefine,'String','0.001');
end
if isfield(Meshopts,'paramaxarea'),
    set(handles.paautomaxarea,'Value',1);
    set(handles.pamaxarea,'Enable','on','String',num2str(Meshopts.paramaxarea));
else
    set(handles.pamaxarea,'String','0.0');
end
if isfield(Meshopts,'primmaxarea'),
    set(handles.prautomaxarea,'Value',1);
    set(handles.prmaxarea,'Enable','on','String',num2str(Meshopts.primmaxarea));
else
    set(handles.prmaxarea,'String','0.001');
end
if isfield(Meshopts,'secmeshrefine'), set(handles.secmeshrefine,'Value',Meshopts.secmeshrefine+1); end

% UIWAIT makes meshopts wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = meshopts_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% --- Executes on button press in pasmooth.
function pasmooth_Callback(hObject, eventdata, handles)


% --- Executes on button press in paautoquality.
function paautoquality_Callback(hObject, eventdata, handles)
wert=get(handles.paautoquality,'Value');
set(handles.paquality,'Enable',onoff(wert));
if wert==0, set(handles.paquality,'String','33.5'); end

% --- Executes on button press in paautorefine.
function paautorefine_Callback(hObject, eventdata, handles)
wert=get(handles.paautorefine,'Value');
set(handles.parefine,'Enable',onoff(wert));
if wert==0, set(handles.parefine,'String','0.2'); end

% --- Executes on button press in paautomaxarea.
function paautomaxarea_Callback(hObject, eventdata, handles)
wert=get(handles.paautomaxarea,'Value');
set(handles.pamaxarea,'Enable',onoff(wert));
if wert==0, set(handles.pamaxarea,'String','0.0'); end

% --- Executes during object creation, after setting all properties.
function paquality_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function paquality_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function parefine_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function parefine_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function pamaxarea_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function pamaxarea_Callback(hObject, eventdata, handles)

% --- Executes on button press in prsmooth.
function prsmooth_Callback(hObject, eventdata, handles)


% --- Executes on button press in prautoquality.
function prautoquality_Callback(hObject, eventdata, handles)
wert=get(handles.prautoquality,'Value');
set(handles.prquality,'Enable',onoff(wert));
if wert==0, set(handles.prquality,'String','33.5'); end

% --- Executes on button press in prautorefine.
function prautorefine_Callback(hObject, eventdata, handles)
wert=get(handles.prautorefine,'Value');
set(handles.prrefine,'Enable',onoff(wert));
if wert==0, set(handles.prrefine,'String','0.001'); end

% --- Executes on button press in prautomaxarea.
function prautomaxarea_Callback(hObject, eventdata, handles)
wert=get(handles.prautomaxarea,'Value');
set(handles.prmaxarea,'Enable',onoff(wert));
if wert==0, set(handles.prmaxarea,'String','0.001'); end

% --- Executes during object creation, after setting all properties.
function prquality_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function prquality_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function prrefine_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function prrefine_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function prmaxarea_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function prmaxarea_Callback(hObject, eventdata, handles)

% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
global Pro Meshopts
ff={'parasmooth','paraquality','pararefine','paramaxarea',...
    'primsmooth','primquality','primrefine','primmaxarea','secmeshrefine'};
for i=1:length(ff), 
    if isfield(Meshopts,ff{i}), Meshopts=rmfield(Meshopts,ff{i}); end 
end
Meshopts.nospline=get(handles.nospline,'Value');
Meshopts.equirefine=get(handles.equirefine,'Value');
if get(handles.pasmooth,'Value')==0, Meshopts.parasmooth=0; end
if get(handles.prsmooth,'Value')==0, Meshopts.primsmooth=0; end
if get(handles.paautoquality,'Value')==1,
    ss=get(handles.paquality,'String');nn=str2num(ss);
    if isnumeric(nn)&&(nn>30)&&(nn<=35), Meshopts.paraquality=nn; end
end
if get(handles.prautoquality,'Value')==1,
    ss=get(handles.prquality,'String');nn=str2num(ss);
    if isnumeric(nn)&&(nn>30)&&(nn<=35), Meshopts.primquality=nn; end
end
if get(handles.paautorefine,'Value')==1,
    ss=get(handles.parefine,'String');nn=str2num(ss);
    if isnumeric(nn), Meshopts.pararefine=nn; end
end
if get(handles.prautorefine,'Value')==1,
    ss=get(handles.prrefine,'String');nn=str2num(ss);
    if isnumeric(nn), Meshopts.primrefine=nn; end
end
if get(handles.paautomaxarea,'Value')==1,
    ss=get(handles.pamaxarea,'String');nn=str2num(ss);
    if isnumeric(nn), Meshopts.paramaxarea=nn; end
end
if get(handles.prautomaxarea,'Value')==1,
    ss=get(handles.prmaxarea,'String');nn=str2num(ss);
    if isnumeric(nn), Meshopts.primmaxarea=nn; end
end
secref=get(handles.secmeshrefine,'Value')-1;
if secref~=1, Meshopts.secmeshrefine=secref; end
delete(gcbf);

% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
delete(gcbf);

% --- Executes on button press in test.
function test_Callback(hObject, eventdata, handles)
global Pro Mesh N progdir
if isfield(Pro,'dirname'), cd(Pro.dirname); end
cmdline=[progdir filesep 'dc2dtreepre -v -t '];
if ~get(handles.nospline,'Value'), cmdline=[cmdline '-B ']; end
if get(handles.equirefine,'Value'), cmdline=[cmdline '-E ']; end
if get(handles.pasmooth,'Value'), cmdline=[cmdline '-S ']; end
if get(handles.paautoquality,'Value')==1, 
    cmdline=[cmdline '-Q' get(handles.paquality,'String') ' ']; end
if get(handles.paautorefine,'Value')==1, 
    cmdline=[cmdline '-R' get(handles.parefine,'String') ' ']; end
if get(handles.paautomaxarea,'Value')==1, 
    cmdline=[cmdline '-A' get(handles.pamaxarea,'String') ' ']; end
cmdline=[cmdline Pro.datfile];
systemcall(cmdline);
% Mesh=loadmesh('tmp/meshPara');
Mesh=loadmesh('tmp\meshPara');
cd(progdir);
Pro.res=ones(Mesh.ncells,1)*median(N.r);
set(figure(1),'MenuBar','none','NumberTitle','off','Name','Parameter Mesh');clf;
tripatchmod(Mesh,Pro.res,N);


function erg=onoff(value)
erg='off';
if value==1, erg='on'; end

% --- Executes during object creation, after setting all properties.
function secmeshrefine_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in secmeshrefine.
function secmeshrefine_Callback(hObject, eventdata, handles)



% --- Executes on button press in nospline.
function nospline_Callback(hObject, eventdata, handles)


% --- Executes on button press in equirefine.
function equirefine_Callback(hObject, eventdata, handles)


