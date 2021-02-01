function varargout = invopts(varargin)
% INVOPTS M-file for invopts.fig
%      INVOPTS, by itself, creates a new INVOPTS or raises the existing
%      singleton*.
%
%      H = INVOPTS returns the handle to a new INVOPTS or the handle to
%      the existing singleton*.
%
%      INVOPTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INVOPTS.M with the given input arguments.
%
%      INVOPTS('Property','Value',...) creates a new INVOPTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before invopts_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to invopts_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help invopts

% Last Modified by GUIDE v2.5 03-Dec-2005 11:50:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @invopts_OpeningFcn, ...
                   'gui_OutputFcn',  @invopts_OutputFcn, ...
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


% --- Executes just before invopts is made visible.
function invopts_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to invopts (see VARARGIN)

% Choose default command line output for invopts
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global Pro Invopts
if isfield(Invopts,'updatesens')&&(Invopts.updatesens>0),
    set(handles.updatesens,'Value',1); end
if isfield(Invopts,'optlambda')&&(Invopts.optlambda>0),
    set(handles.optlambda,'Value',1); end
if isfield(Invopts,'robustdata')&&(Invopts.robustdata>0),
    set(handles.robustdata,'Value',1); end
if isfield(Invopts,'blockymodel')&&(Invopts.blockymodel>0),
    set(handles.blockymodel,'Value',1); end
if isfield(Invopts,'lowerbound')&&isnumeric(Invopts.lowerbound)&&(Invopts.lowerbound>0),
    set(handles.lowerbound,'String',num2str(Invopts.lowerbound)); end
if isfield(Invopts,'upperbound')&&isnumeric(Invopts.upperbound)&&(Invopts.upperbound>0),
    set(handles.upperbound,'String',num2str(Invopts.upperbound)); end
% UIWAIT makes invopts wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = invopts_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% --- Executes on button press in ok.
function ok_Callback(hObject, eventdata, handles)
global Invopts
Invopts=[];
cc=get(handles.constraints,'Value')-1;
if cc>0, Invopts.constraints=cc; end
ll=str2num(get(handles.lowerbound,'String'));
if isnumeric(ll)&&(ll>0), Invopts.lowerbound=ll; else ll=0; end
cc=str2num(get(handles.upperbound,'String'));
if isnumeric(cc)&&(cc>ll), Invopts.upperbound=cc; end
cc=get(handles.optlambda,'Value');
if cc, Invopts.optlambda=1; end
cc=get(handles.updatesens,'Value');
if cc, Invopts.updatesens=1; end
cc=get(handles.blockymodel,'Value');
if cc, Invopts.blockymodel=1; end
cc=get(handles.robustdata,'Value');
if cc, Invopts.robustdata=1; end
delete(gcbf);

% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
delete(gcbf);

% --- Executes on button press in test.
function test_Callback(hObject, eventdata, handles)
global Pro S N
if isfield(Pro,'dirname'), cd(Pro.dirname); end
% systemcall(cmdline);
% % Mesh=loadmesh('tmp/meshPara');
% Mesh=loadmesh('tmp\meshPara');
% cd(progdir);
% Pro.res=ones(Mesh.ncells,1)*median(N.r);
% figure(1);clf;
% tripatchmod(Mesh,Pro.res,N);


function erg=onoff(value)
erg='off';
if value==1, erg='on'; end

% --- Executes during object creation, after setting all properties.
function constraints_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in constraints.
function constraints_Callback(hObject, eventdata, handles)


% --- Executes on button press in updatesens.
function updatesens_Callback(hObject, eventdata, handles)


% --- Executes on button press in optlambda.
function optlambda_Callback(hObject, eventdata, handles)


% --- Executes on button press in robustdata.
function robustdata_Callback(hObject, eventdata, handles)


% --- Executes on button press in blockymodel.
function blockymodel_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function upperbound_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function upperbound_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function lowerbound_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function lowerbound_Callback(hObject, eventdata, handles)

