function varargout = s_par(varargin)
% S_PAR M-file for s_par.fig
%      S_PAR, by itself, creates a new S_PAR or raises the existing
%      singleton*.
%
%      H = S_PAR returns the handle to a new S_PAR or the handle to
%      the existing singleton*.
%
%      S_PAR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in S_PAR.M with the given input arguments.
%
%      S_PAR('Property','Value',...) creates a new S_PAR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before s_par_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to s_par_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help s_par

% Last Modified by GUIDE v2.5 15-Feb-2004 16:09:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @s_par_OpeningFcn, ...
                   'gui_OutputFcn',  @s_par_OutputFcn, ...
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


% --- Executes just before s_par is made visible.
function s_par_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for s_par
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
global Mod P
lxy=(length(Mod.x)-1)*(length(Mod.z)-1);
if size(P,1)~=lxy, P=speye(lxy); end
set(handles.xa,'String',num2str(min(Mod.x)));
set(handles.xe,'String',num2str(max(Mod.x)));
set(handles.za,'String',num2str(min(Mod.z)));
set(handles.ze,'String',num2str(max(Mod.z)));
set(handles.from,'String',num2str(size(P,1)));
set(handles.to,'String',num2str(size(P,2)));
% UIWAIT makes s_par wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = s_par_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in combine.
function combine_Callback(hObject, eventdata, handles)

% --- Executes on button press in delete.
function delete_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function xa_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function xa_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function xe_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function xe_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function edit3_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function ze_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function ze_Callback(hObject, eventdata, handles)

% --- Executes on button press in add.
function add_Callback(hObject, eventdata, handles)
global Mod
xq=Mod.x(1:end-1)+diff(Mod.x)/2;zq=Mod.z(1:end-1)+diff(Mod.z)/2;
xa=str2num(get(handles.xa,'String'));
xe=str2num(get(handles.xe,'String'));
za=str2num(get(handles.za,'String'));
ze=str2num(get(handles.ze,'String'));
ia=max(find([-Inf Mod.x]<=xa))-1;ie=min(find([Mod.x Inf]>=xe))-1;
ka=max(find([-Inf Mod.z]<=za))-1;ke=min(find([Mod.z Inf]>=ze))-1;
if (xa==xe)&&any(Mod.x==xa), ie=ie+1;ia=ia-1; end
if (za==ze)&&any(Mod.z==za), ke=ke+1;ka=ka-1; end
% [ia ie ka ke]
% [Mod.x(ia) Mod.x(ie+1) Mod.z(ka) Mod.z(ke+1)]
% return
ii=ie-ia+1;kk=ke-ka+1;
if ii*kk>0,
    ss=get(handles.regions,'String');
    ss{end+1}=sprintf('%4g %4g %4g %4g %dx%d=%d cells',...
        xa,xe,za,ze,ii,kk,ii*kk);
    set(handles.regions,'String',ss);
end

% --- Executes during object creation, after setting all properties.
function za_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function za_Callback(hObject, eventdata, handles)

% --- Executes on button press in apply.
function apply_Callback(hObject, eventdata, handles)
global Mod P Cov INV
mincov=0.4;if isfield(INV,'mico'), mincov=INV.mico; end
lxy=(length(Mod.x)-1)*(length(Mod.z)-1);lx=length(Mod.x)-1;
if get(handles.combine,'Value'),
    P=pmatrix2d(Mod.x,Mod.z);
else
    P=speye(lxy);
end
if get(handles.delete,'Value'),
    pcov=Cov(:)'*P./sum(P);
    P(:,find(pcov<mincov))=[];
end
ss=get(handles.regions,'String');
for i=3:length(ss),
    aa=sscanf(ss{i},'%g %g %g %g');
    xa=aa(1);xe=aa(2);za=aa(3);ze=aa(4);
    ia=max(find([-Inf Mod.x]<=xa))-1;ie=min(find([Mod.x Inf]>=xe))-1;
    ka=max(find([-Inf Mod.z]<=za))-1;ke=min(find([Mod.z Inf]>=ze))-1;
    if (xa==xe)&find(Mod.x==xa), ie=ie+1;ia=ia-1; end
    if (za==ze)&find(Mod.z==za), ke=ke+1;ka=ka-1; end
    fi=[];
    for k=ka:ke,
        for ii=ia:ie,
            ind=(k-1)*lx+ii;
            fi=[fi find(P(ind,:))];
        end
    end
    P(:,fi)=[];
end
fpar=size(P,2);
set(handles.to,'String',num2str(fpar));


% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
if ishandle(99), delete(99); end
delete(gcbf);

% --- Executes on button press in show.
function show_Callback(hObject, eventdata, handles)
global Mod P
mm=zeros(length(Mod.x)-1,length(Mod.z)-1);
fpar=size(P,2);
mm(:)=0;mm(:)=P*(rand(fpar,1)*10+10);
figure(99);draw2dmodel(Mod.x,Mod.z,mm);
cmap=colormap;cmap(1,:)=1;colormap(cmap);

