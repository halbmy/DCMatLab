function varargout = daterr(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @daterr_OpeningFcn, ...
                   'gui_OutputFcn',  @daterr_OutputFcn, ...
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


% --- Executes just before daterr is made visible.
function daterr_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for daterr
handles.output = hObject;
guidata(hObject, handles);
iconify(hObject);
global N
if isfield(N,'i')&&(length(N.i)==length(N.r)), 
    set(handles.current,'String','meas.'); 
    set(handles.current,'Enable','Off'); 
end
if isfield(N,'proz'), set(handles.proz,'String',num2str(N.proz)); end
if isfield(N,'umin'), set(handles.umin,'String',num2str(N.umin)); end
if isfield(N,'i')&&(length(N.i)==1), 
    set(handles.current,'String',num2str(N.i)); 
end
daterr('update_Callback',hObject,[],guidata(hObject))
if ~isfield(N,'err')||isempty(N.err),
    set(handles.add,'Enable','off');
    set(handles.ok,'Enable','off');
    set(handles.show,'Enable','off');
end

% UIWAIT makes daterr wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = daterr_OutputFcn(hObject, eventdata, handles)
%varargout{1} = handles.output;

% --- Executes on button press in const.
function const_Callback(hObject, eventdata, handles)
global N
N.err=ones(size(N.r))*0.03;
uiresume(handles.figure1);
delete(gcbf);

% --- Executes on button press in readfromfile.
function readfromfile_Callback(hObject, eventdata, handles)
global N
errfile=uigetfile('*.*','Read error file');
if exist(errfile,'file'),
    N.err=load(errfile,'-ascii');
else
    return
end
set(handles.status,'String',sprintf('min=%.1f%% max=%.1f%%',min(N.err)*100,max(N.err)*100));
% pause(1.0);
uiresume(handles.figure1);
delete(gcbf);

% --- Executes during object creation, after setting all properties.
function proz_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function proz_Callback(hObject, eventdata, handles)
daterr('update_Callback',gcbo,[],guidata(gcbo))

% --- Executes during object creation, after setting all properties.
function umin_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function umin_Callback(hObject, eventdata, handles)
daterr('update_Callback',gcbo,[],guidata(gcbo))

% --- Executes during object creation, after setting all properties.
function current_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function current_Callback(hObject, eventdata, handles)
daterr('update_Callback',gcbo,[],guidata(gcbo))

function update_Callback(hObject, eventdata, handles)
global N
proz=str2double(get(handles.proz,'String'));
umin=str2double(get(handles.umin,'String'));
if isfield(N,'i')&&(length(N.i)==length(N.r)),
    current=N.i;
else
    current=str2double(get(handles.current,'String'))/1000;
end
err=estimateerror(N,proz,umin/1000,current);
elvar=str2double(get(handles.elvar,'String')); 
if elvar>0,
    [A,maxerr]=electrodeerr3d(N,elvar);
    err=err+0.5*maxerr;
end 
set(handles.status,'String',sprintf('min=%.1f%% max=%.1f%% mean=%.1f%%',min(err)*100,max(err)*100,mean(err)*100));
handles.output=err;

% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
uiresume(handles.figure1);
delete(gcbf);

% --- Executes on button press in show.
function show_Callback(hObject, eventdata, handles)
global N MAL
% proz=str2double(get(handles.proz,'String'));
% umin=str2double(get(handles.umin,'String'));
% if isfield(N,'i')&&(length(N.i)==length(N.r)),
%     current=N.i;
% else
%     current=str2double(get(handles.current,'String'))/1000;
% end
% err=estimateerror(N,proz,umin/1000,current);
% elvar=str2double(get(handles.elvar,'String')); 
% if elvar>0,
%     maxerr=electrodeerr3d(N,elvar);
%     err=err+maxerr;
% end 
mal=MAL;mal.cauto=1;
set(figure(2),'NumberTitle','off','Name','Error');
iconify(2);
mal=MAL;mal.vie=1;mal.clog=1;mal.cauto=1;
[cmin,cmax]=plotprofiles(N,N.err*100,mal);

% --- Executes during object creation, after setting all properties.
function elvar_CreateFcn(hObject, eventdata, handles)

function elvar_Callback(hObject, eventdata, handles)
daterr('update_Callback',gcbo,[],guidata(gcbo))


% --- Executes on button press in set.
function set_Callback(hObject, eventdata, handles)
global N
N.proz=str2double(get(handles.proz,'String'));
N.umin=str2double(get(handles.umin,'String'));
if isfield(N,'i')&&(length(N.i)==length(N.r)),
    current=N.i;
else
    current=str2double(get(handles.current,'String'))/1000;
    N.i=current;
end
N.err=estimateerror(N,N.proz,N.umin/1000,current);
elvar=str2double(get(handles.elvar,'String')); 
if elvar>0,
    [A,maxerr]=electrodeerr3d(N,elvar);
    N.err=N.err+0.5*maxerr;
end 
set(handles.add,'Enable','on');
set(handles.ok,'Enable','on');
set(handles.show,'Enable','on');

% --- Executes on button press in add.
function add_Callback(hObject, eventdata, handles)
global N
proz=str2double(get(handles.proz,'String'));
umin=str2double(get(handles.umin,'String'));
if isfield(N,'i')&&(length(N.i)==length(N.r)),
    current=N.i;
else
    current=str2double(get(handles.current,'String'))/1000;
end
N.err=N.err+estimateerror(N,proz,umin/1000,current);
elvar=str2double(get(handles.elvar,'String')); 
if elvar>0,
    [A,maxerr]=electrodeerr3d(N,elvar);
    N.err=N.err+0.5*maxerr;
end 
set(handles.status,'String',sprintf('min=%.1f%% max=%.1f%%',min(N.err)*100,max(N.err)*100));

