function varargout = debyegui(varargin)
% DEBYEGUI M-file for debyegui.fig
%      DEBYEGUI, by itself, creates a new DEBYEGUI or raises the existing
%      singleton*.
%
%      H = DEBYEGUI returns the handle to a new DEBYEGUI or the handle to
%      the existing singleton*.
%
%      DEBYEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEBYEGUI.M with the given input arguments.
%
%      DEBYEGUI('Property','Value',...) creates a new DEBYEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before debyegui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to debyegui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help debyegui

% Last Modified by GUIDE v2.5 21-Oct-2010 14:32:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @debyegui_OpeningFcn, ...
                   'gui_OutputFcn',  @debyegui_OutputFcn, ...
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


% --- Executes just before debyegui is made visible.
function debyegui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to debyegui (see VARARGIN)

% Choose default command line output for debyegui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes debyegui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = debyegui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function File_Callback(hObject, eventdata, handles)
% hObject    handle to File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function FileOpen_Callback(hObject, eventdata, handles)
alltypes='*.res;*.RES;*.dat';infile=alltypes;
if isfield(handles,'datafile'), 
%     handles.datafile
    infile=fullfile(fileparts(handles.datafile),infile);
end
[fname,pname]=uigetfile(alltypes,'Read data file',infile);
if ~ischar(fname), return; end
handles.datafile=fullfile(pname,fname);
Data=read1resfile(handles.datafile);
if ~isstruct(Data)||~isfield(Data,'f')||isempty(Data.f), return; end
handles.data=Data;
guidata(hObject,handles);
taumax=rndig(1/min(Data.f)/2/pi*8);
set(handles.taumax,'String',num2str(taumax));
set(handles.fmin,'String','0');
axes(handles.axes2);cla;
debyegui('ShowRawData_Callback',hObject,eventdata,guidata(hObject))

% --------------------------------------------------------------------
function FileExit_Callback(hObject, eventdata, handles)
delete(gcbo);


% --- Executes on button press in openfilebutton.
function openfilebutton_Callback(hObject, eventdata, handles)
debyegui('FileOpen_Callback',hObject,eventdata,guidata(hObject))

% --------------------------------------------------------------------
function Show_Callback(hObject, eventdata, handles)
% hObject    handle to Show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function ShowRawData_Callback(hObject, eventdata, handles)
f=handles.data.f;
phi=handles.data.phi;
dphi=handles.data.dphi;
kfak=str2num('1/39.8');
rho=handles.data.rhoa*kfak;
drho=handles.data.drhoa*kfak;
axes(handles.axes2);
semilogx(f,rho,'bx-');
grid on;
for i=1:length(f), line([1 1]*f(i),([-1 1]*drho(i)/100+1)*rho(i),'Color','blue'); end
xlim([min(f)*0.9 max(f)*1.1]);
xlabel('f in Hz');
ylabel('\rho in \Omegam');
axes(handles.axes1);
% errorbar(log10(f),phi,phi-dphi,phi+dphi);
if min(phi)>0,
    loglog(f,phi*1000,'bx-');
else
    semilogx(f,phi*1000,'bx-');
end
for i=1:length(f), line([1 1]*f(i),[-1 1]*dphi(i)*1000+phi(i)*1000,'Color','blue'); end
grid on;
xlim([min(f)*0.9 max(f)*1.1]);
xlabel('f in Hz');
ylabel('\phi in mrad');


% --- Executes on button press in removeembutton.
function removeembutton_Callback(hObject, eventdata, handles)
lam=str2num(get(handles.lamem,'String'));
fdamp=str2num(get(handles.fdamp,'String'));
opt=struct('lam',lam,'fdamp',fdamp);
f=handles.data.f;
phi=handles.data.phi;
dphi=handles.data.dphi;
phi1=removeemcoupling(f,phi,dphi,opt);
axes(handles.axes2);
semilogx(f,phi1*1000,'bx-');
for i=1:length(f), line([1 1]*f(i),[-1 1]*dphi(i)*1000+phi1(i)*1000,'Color','blue'); end
grid on;
xlim([min(f)*0.9 max(f)*1.1]);
xlabel('f in Hz');
ylabel('\phi in mrad');
handles.data.phi1=phi1;
guidata(hObject,handles);

function fdamp_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function fdamp_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lamem_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of lamem as text
%        str2double(get(hObject,'String')) returns contents of lamem as a double


% --- Executes during object creation, after setting all properties.
function lamem_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function ShowData_Callback(hObject, eventdata, handles)
axes(handles.axes1);
f=handles.data.f;
phi=handles.data.phi;
if isfield(handles.data,'phi1'), phi=handles.data.phi1; end
dphi=handles.data.dphi;
% errorbar(log10(f),phi,phi-dphi,phi+dphi);
semilogx(f,phi*1000,'bx-');
for i=1:length(f), line([1 1]*f(i),[-1 1]*dphi(i)*1000+phi(i)*1000,'Color','blue'); end
grid on;
xlim([min(f)*0.9 max(f)*1.1]);
xlabel('f in Hz');
ylabel('\phi in mrad');

% --------------------------------------------------------------------
function ShowDebye_Callback(hObject, eventdata, handles)
axes(handles.axes2);cla;
tau=handles.model.tau;
m=handles.model.m;
semilogx(tau,m,'bx-');
grid on;
xlim(minmax(tau));
xlabel('\tau in s');
ylabel('spectral chargeability');
yl=ylim;
t50=exp(sum(log(tau).*m)/sum(m));
line(t50*[1 1],yl,'Color','r');
t=text(t50,yl(2),num2str(rndig(t50*1000),'%g ms'),'Color','red');
set(t,'HorizontalAlignment','center','VerticalAlignment','bottom');
[ma,loc]=max(m);tmax=tau(loc);
line(tmax*[1 1],yl,'Color','b');
t=text(tmax,yl(2),num2str(rndig(tmax*1000),'%g ms'),'Color','blue');
set(t,'HorizontalAlignment','center','VerticalAlignment','bottom');



function fmax_Callback(hObject, eventdata, handles)
fmax=str2num(get(hObject,'String'));
taumin=rndig(1/fmax/2/pi/8);
set(handles.taumin,'String',num2str(taumin));

% --- Executes during object creation, after setting all properties.
function fmax_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function taumin_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of taumin as text
%        str2double(get(hObject,'String')) returns contents of taumin as a double


% --- Executes during object creation, after setting all properties.
function taumin_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function taumax_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of taumax as text
%        str2double(get(hObject,'String')) returns contents of taumax as a double


% --- Executes during object creation, after setting all properties.
function taumax_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ntau_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of ntau as text
%        str2double(get(hObject,'String')) returns contents of ntau as a double


% --- Executes during object creation, after setting all properties.
function ntau_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lam_Callback(hObject, eventdata, handles)
debyegui('dodebye_Callback',hObject,eventdata,guidata(hObject))


% --- Executes during object creation, after setting all properties.
function lam_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dodebye.
function dodebye_Callback(hObject, eventdata, handles)
f=handles.data.f;
phi=handles.data.phi;
if isfield(handles.data,'phi1'), phi=handles.data.phi1; end
dphi=handles.data.dphi;
fmin=str2num(get(handles.fmin,'String'));
fmax=str2num(get(handles.fmax,'String'));
fi=find((f<=fmax)&(f>=fmin));
f=f(fi);
phi=phi(fi);
dphi=dphi(fi);
lam=str2num(get(handles.lam,'String'));
ntau=str2num(get(handles.ntau,'String'));
taumin=str2num(get(handles.taumin,'String'));
taumax=str2num(get(handles.taumax,'String'));
tau=logspace(log10(taumin),log10(taumax),ntau)'; % tau discretisation
[m,phiout]=debyedecomp(f,phi,dphi,tau,lam,-1);
handles.model.m=m;
handles.model.tau=tau;
guidata(hObject,handles);
debyegui('ShowData_Callback',hObject,eventdata,guidata(hObject))
hold on;
plot(f,phiout*1000,'ro-');
hold off;
debyegui('ShowDebye_Callback',hObject,eventdata,guidata(hObject))


% --- Executes on button press in lamup.
function lamup_Callback(hObject, eventdata, handles)
slam=get(handles.lam,'String');
lam=rndig(str2num(slam)*1.5);
set(handles.lam,'String',num2str(lam));
debyegui('lam_Callback',hObject,eventdata,guidata(hObject))


% --- Executes on button press in lamdown.
function lamdown_Callback(hObject, eventdata, handles)
slam=get(handles.lam,'String');
lam=rndig(str2num(slam)/1.5);
set(handles.lam,'String',num2str(lam));
debyegui('lam_Callback',hObject,eventdata,guidata(hObject))


% --------------------------------------------------------------------
function Export_Callback(hObject, eventdata, handles)
isex='Off';
if isfield(handles,'datafile'), isex='on'; end
set(handles.ExportFigure,'Enable',isex);
if ~isfield(handles,'data'), isex='Off'; end
set(handles.ExportSpectrum,'Enable',isex);
set(handles.ExportData,'Enable',isex);
if ~isfield(handles,'model'), isexe='off'; end
set(handles.ExportChargeability,'Enable',isex);
set(handles.ExportModel,'Enable',isex);

% --------------------------------------------------------------------
function ExportFigure_Callback(hObject, eventdata, handles)
[pp,nn,ee]=fileparts(handles.datafile);
outfile=strrep(handles.datafile,ee,'.fig');
[fname,fpath]=uiputfile('*.fig','Export Figure File',outfile);
if isstr(fname),
   outfile=fullfile(fpath,fname);
   saveas(handles.figure1,outfile,'fig');
end

% --------------------------------------------------------------------
function ExportSpectrum_Callback(hObject, eventdata, handles)
figure(11);set(11,'Units','characters','Position',get(handles.figure1,'Position'));
copyobj(handles.axes1,11);
[pp,nn,ee]=fileparts(handles.datafile);
outfile=strrep(handles.datafile,ee,'-data.pdf');
[fname,fpath]=uiputfile('*.pdf','Export Data Plot',outfile);
if isstr(fname),
   outfile=fullfile(fpath,fname);
   epsprint(11,strrep(outfile,'.pdf','.eps'),1);
end
close(11);


% --------------------------------------------------------------------
function ExportChargeability_Callback(hObject, eventdata, handles)
figure(11);set(11,'Units','characters','Position',get(handles.figure1,'Position'));
copyobj(handles.axes2,11);
[pp,nn,ee]=fileparts(handles.datafile);
outfile=strrep(handles.datafile,ee,'-model.pdf');
[fname,fpath]=uiputfile('*.pdf','Export Model Plot',outfile);
if isstr(fname),
   outfile=fullfile(fpath,fname);
   epsprint(11,strrep(outfile,'.pdf','.eps'),1);
end
close(11);


% --------------------------------------------------------------------
function ExportData_Callback(hObject, eventdata, handles)
[pp,nn,ee]=fileparts(handles.datafile);
outfile=strrep(handles.datafile,ee,'-data.txt');
[fname,fpath]=uiputfile('*.txt','Data Text File',outfile);
if isstr(fname),
   outfile=fullfile(fpath,fname);
   fid=fopen(outfile,'w');
   fprintf(fid,'%8s\t%8s\t%8s\t%8s\t%8s\r\n','f/Hz','Ra/Ohmm','phi/mrad','dR/%','dphi/mrad');
   Data=handles.data;phi=Data.phi;
   if isfield(Data,'phi1'), phi=Data.phi1; end
   A=[Data.f(:) Data.rhoa(:) phi(:)*1000 Data.drhoa(:) Data.dphi(:)*1000];
   fprintf(fid,'%8g\t%8g\t%8g\t%8g\t%8g\r\n',A');
   fclose(fid);
end


% --------------------------------------------------------------------
function ExportModel_Callback(hObject, eventdata, handles)
[pp,nn,ee]=fileparts(handles.datafile);
outfile=strrep(handles.datafile,ee,'-model.txt');
[fname,fpath]=uiputfile('*.txt','Model Text File',outfile);
if isstr(fname),
   outfile=fullfile(fpath,fname);
   fid=fopen(outfile,'w');
   fprintf(fid,'%8s\t%8s\r\n','tau/s','m/1e-6');
   A=[handles.model.tau(:) handles.model.m(:)*1e6];
   fprintf(fid,'%8g\t%8g\r\n',A');
   fclose(fid);
end



function fmin_Callback(hObject, eventdata, handles)
fmin=str2num(get(hObject,'String'));
if fmin>0,
    taumax=rndig(1/fmin/pi);
    % set(handles.taumax,'String',num2str(taumax));
end


% --- Executes during object creation, after setting all properties.
function fmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
