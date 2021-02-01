function varargout = tt2dtomo(varargin)
% TT2DTOMO M-file for tt2dtomo.fig
%      TT2DTOMO, by itself, creates a new TT2DTOMO or raises the existing
%      singleton*.
%
%      H = TT2DTOMO returns the handle to a new TT2DTOMO or the handle to
%      the existing singleton*.
%
%      TT2DTOMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TT2DTOMO.M with the given input arguments.
%
%      TT2DTOMO('Property','Value',...) creates a new TT2DTOMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tt2dtomo_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tt2dtomo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tt2dtomo

% Last Modified by GUIDE v2.5 23-Mar-2006 11:08:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tt2dtomo_OpeningFcn, ...
                   'gui_OutputFcn',  @tt2dtomo_OutputFcn, ...
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


% --- Executes just before tt2dtomo is made visible.
function tt2dtomo_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for tt2dtomo
handles.output = hObject;
set(handles.text14,'String','TT2dTomo v0.3')
set(handles.figure1,'Name','TT2dTomo - Refraction tomography by resistivity.net');
% Update handles structure
guidata(hObject, handles);
% set(handles.upper,'Enable','on');
% set(handles.dx,'String',0.5);
% UIWAIT makes tt2dtomo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tt2dtomo_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in file.
function file_Callback(hObject, eventdata, handles)
global Shot datfile err
% [fname,pname]=uigetfile('*.txt;*.grm;*.tom','Data File');
% if ~ischar(fname), return; end
allfiles='*.pit;*.txt;*.tom;*.dat';
if isempty(datfile),
    infile=allfiles;
else
    [pp,nn,ee]=fileparts(datfile);
    infile=fullfile(pp,allfiles);
end
[fname,pname]=uigetfile(allfiles,'Data File',infile);
if isequal(fname,0), return; end
datfile=fullfile(pname,fname);
[pp,ff,ee]=fileparts(datfile);
if isequal(lower(ee),'.pit'),
    Shot=readpit(datfile);
elseif isequal(lower(ee),'.dat'),
    fid=fopen(datfile,'r');zeile=fgetl(fid);fclose(fid);
    if length(str2num(zeile))>3, 
%         Shot=korrrad(readjuett(datfile));
        Shot=readjuett(datfile);
        Shot.pos=round(Shot.pos*1000)/1000;
    else Shot=readunishot(datfile); end
elseif isequal(lower(ee),'.tom'),
    Shot=readtom(datfile);
    Shot.t=Shot.t/1000;
else
    errordlg('File type not recognized!');return;
end
set(handles.message,'String',...
    sprintf('%d data with %d points',length(Shot.t),size(Shot.pos,1)));
if ~isfield(Shot,'sd'), Shot.sd=0.05*Shot.t; end
if ~isfield(Shot,'va'),
    Shot.va=[];l=0;
    for i=1:length(Shot.ns),
       aa=Shot.pos(Shot.nx{i},:);
       aa(:,1)=aa(:,1)-Shot.pos(Shot.ns{i},1);
       aa(:,2)=aa(:,2)-Shot.pos(Shot.ns{i},2);
       dist=sqrt(sum(aa.^2,2));
       Shot.va=[Shot.va;dist./Shot.t(l+(1:length(dist)))];
       l=l+length(dist);
    end
end
err=Shot.sd;
if max(Shot.t)>1e-3, err(err<1e-5)=1e-5; end
axes(handles.malfeld);
shotimage(Shot);
nn=[ 'TT2DTomo - ' strrep(datfile,[pwd filesep],'')];
try, set(hObject,'Name',nn);
catch, set(handles.figure1,'Name',nn); end

% --- Executes on button press in makemesh.
function makemesh_Callback(hObject, eventdata, handles)
global Shot Mesh datfile velocity W t chiq C err
set(gcf,'Pointer','Watch');
err=Shot.sd;
if max(Shot.t)>1e-3, err(err<1e-5)=1e-5; end
[pp,ff,ee]=fileparts(datfile);
meshname=strrep(datfile,ee,'');
dx=str2num(get(handles.dx,'String'));
if ~isnumeric(dx), dx=0.1; end
if 0,
    writeroundpoly(Shot.pos,'test.d',dx,1);
    dos('easymesh test.d');
    Mesh=readeasymesh('test');
else
    writeroundpoly(Shot.pos,[meshname '.poly'],dx);
    maxarea=sum(diff(Shot.pos(1:2,:)).^2)*dx^2*4;
    if maxarea==0, maxarea=0.003; end
    if dx==0, maxarea=0; end
    quality=33.8;
    if maxarea>0,
        dos(['dctriangle -S -q' num2str(quality) ' -a' num2str(maxarea) ' "' meshname '.poly"']);
    else
        dos(['dctriangle -S -q' num2str(quality) ' "' meshname '.poly"']);
    end
    Mesh=loadmesh([meshname '.bms']);
end
axes(handles.malfeld);cla reset;
patch('Vertices',Mesh.node,'Faces',Mesh.cell,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]);
axis equal tight
nsm=cell2mat(Shot.ns);hold on;plot(Shot.pos(nsm,1),Shot.pos(nsm,2),'k.');hold off
drawnow;
v0=median(Shot.va);
velocity=ones(Mesh.ncells,1)*v0;
[W,t]=waymatrix(Mesh,Shot,velocity);
C=onesidesmoothness(Mesh,0);
chiq=chi2(Shot.t,t,err,0);
% set(handles.message,'String',...
%     sprintf('CHI^2=%.2f SQD=%.2fµs\n',chiq(end),sqrt(mean((Shot.t-t).^2))*1e6));
set(handles.message,'String',...
    sprintf('CHI^2=%.2f SQD=%ss\n',chiq(end),int2hum(sqrt(mean((Shot.t-t).^2)),1)));
mal=struct('cauto',0,'cmin',v0/1.1,'cmax',v0*1.1);
cla reset;tripatchmod(Mesh,velocity,1,mal);
set(gcf,'Pointer','Arrow');

function lower_Callback(hObject, eventdata, handles)
if ~isnumeric(str2num(get(hObject,'String'))), set(hObject,'String','0'); end

% --- Executes during object creation, after setting all properties.
function lower_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function lambda_Callback(hObject, eventdata, handles)
if ~isnumeric(str2num(get(hObject,'String'))), set(hObject,'String','0'); end

% --- Executes during object creation, after setting all properties.
function lambda_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
global Mesh velocity Shot t W C chiq err
deltaData=Shot.t-t;
D=spdiags(1./err,0,length(err),length(err));
D2=D*D;lv=length(velocity);
lbound=str2num(get(handles.lower,'String'));
if lbound>min(velocity), lbound=0;set(handles.lower,'String','0'); end
ubound=0;
lam=str2double(get(handles.lambda,'String'));
set(gcf,'Pointer','Watch');
deltaModel=log(velocity-lbound);
vv=-(velocity-lbound)./velocity.^2;
S=W*spdiags(vv,0,lv,lv);
% ss=svd(full(S));max(ss)
%     dm=cglscdp(S,deltaData,lam,C,D,1,deltaModel);
dm=(S'*D2*S+lam*C'*C)\(S'*(D2*deltaData)-lam*(C'*C*deltaModel));
if ubound>0,
    velocity1=(ubound*(velocity-lbound).*exp(dm)+lbound*(ubound-velocity))./((velocity-lbound).*exp(dm)+ubound-velocity);
else
    velocity1=(velocity-lbound).*exp(dm)+lbound;
end
fprintf('min/max v=%g/%g',rndig(min(velocity1)),rndig(max(velocity1)));
mal=struct('canot','m/s','clog',1);
axes(handles.malfeld);
tripatchmod(Mesh,velocity1,1,mal);drawnow;
oldt=t;[W,t]=waymatrix(Mesh,Shot,velocity1);
for i=1:21, % line search procedure
    tau=0.05*i;ta=tau*t+(1-tau)*oldt;
    phia(i)=chi2(Shot.t,ta,err,0);
end
[pp,ii]=min(phia);fak=ii*0.05;
if fak<1,
    fprintf('  Line search factor %.2f',fak);dm=dm*fak;
    if ubound>0,
        velocity=(ubound*(velocity-lbound).*exp(dm)+lbound*(ubound-velocity))./((velocity-lbound).*exp(dm)+ubound-velocity);
    else
        velocity=(velocity-lbound).*exp(dm)+lbound;
    end
    fprintf('  min/max v=%d/%d',round(min(velocity)),round(max(velocity)));
    [W,t]=waymatrix(Mesh,Shot,velocity);
else
    velocity=velocity1;
end
fprintf('\n');
% coverage=sum(W);scoverage=(C'*(C*coverage(:))~=0);
scoverage=1;
axes(handles.malfeld);
tripatchmod(Mesh,velocity,scoverage,mal);%axis normal
axes(handles.malfeld);
nsm=cell2mat(Shot.ns);hold on;
plot(Shot.pos(nsm,1),Shot.pos(nsm,2),'k.');hold off
set(gca,'xlim',[min(Mesh.node(:,1)) max(Mesh.node(:,1))]);
chiq(end+1)=chi2(Shot.t,t,err,0);
% set(handles.message,'String',...
%     sprintf('CHI^2=%.2f SQD=%.2fµs\n',chiq(end),sqrt(mean((Shot.t-t).^2))*1e6));
set(handles.message,'String',...
    sprintf('CHI^2=%.2f SQD=%ss\n',chiq(end),int2hum(sqrt(mean((Shot.t-t).^2)))));
if get(handles.blocky,'Value'), % blocky inversion
    dx=str2num(get(handles.dx,'String'));
    C1=onesidesmoothness(Mesh,dx);
    sm=abs(C1*log(velocity));su2=sum(sm.^2);sua=sum(sm);
    wxz=ones(size(sm));fi=find(sm);
    if ~isempty(fi), wxz(fi)=su2/sua./sm(fi); end
    wxz(wxz>1)=1;
    C=spdiags(wxz,0,length(wxz),length(wxz))*C1;
end
if get(handles.robust,'Value')>0,
    dt=(Shot.t-t)./err; % robust data inversion
    w=abs(dt)*sum(abs(dt))/sum(dt.^2);
    w(w<1)=1;err=err.*w;
end
set(gcf,'Pointer','Arrow');


% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
delete(gcbf);


% --- Executes during object creation, after setting all properties.
function zpower_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function zpower_Callback(hObject, eventdata, handles)
global Mesh C
if ~isnumeric(str2num(get(hObject,'String'))), set(hObject,'String','0'); end
zpower=str2num(get(handles.zpower,'String'));
C=onesidesmoothness(Mesh,zpower);

% --- Executes during object creation, after setting all properties.
function dx_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function dx_Callback(hObject, eventdata, handles)
dx=str2num(get(hObject,'String'));
if ~isnumeric(dx), set(hObject,'String',0); end

% --- Executes on button press in blocky.
function blocky_Callback(hObject, eventdata, handles)

% --- Executes on button press in robust.
function robust_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function fileopen_Callback(hObject, eventdata, handles)
tt2dtomo('file_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function loadmesh_Callback(hObject, eventdata, handles)
global datfile Mesh
[pp,nn,ee]=fileparts(datfile);
[fname,pname]=uigetfile({'*.bms','Mesh files'},'Select mesh File',[pp filesep '*.bms']);
Mesh=loadmesh(fullfile(pname,fname));

% --------------------------------------------------------------------
function fileexit_Callback(hObject, eventdata, handles)
delete(gcbf);

% --------------------------------------------------------------------
function options_Callback(hObject, eventdata, handles)
set(handles.optmesh,'Enable','Off');
set(handles.optgraphics,'Enable','Off');
set(handles.optinv,'Enable','Off');

% --------------------------------------------------------------------
function optmesh_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function optgraphics_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function optinv_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function exportmenu_Callback(hObject, eventdata, handles)
set(handles.exmodel,'Enable','off');
set(handles.exdata,'Enable','off');

% --------------------------------------------------------------------
function exfigure_Callback(hObject, eventdata, handles)
global datfile
[pp,nn,ee]=fileparts(datfile);
[fname,pname]=uiputfile('*.eps','Export Graphics',strrep(datfile,ee,'.eps'));
if ischar(fname),
%     print(handles.figure1,'-depsc2',fullfile(pname,fname));
    epsprint(handles.figure1,fullfile(pname,fname));
end


% --------------------------------------------------------------------
function exmodel_Callback(hObject, eventdata, handles)
% set(handles.malfeld,'Units','pixel');
% poall=get(handles.inv2d,'Position');
% rel=get(handles.message,'Position');rel=rel(4)/poall(4); %new!
% pof=get(gcbf,'position');pof(4)=pof(4)*(1-rel);
% set(f,'Units','Character','position',pof); % gleiche figure-groesse
% set(handles.malfeld,'Units','normalized'); % zurück
% po=get(handles.malfeld,'position');
% set(f,'Units','pixel','PaperSize',po(3:4),'PaperPositionMode','auto');
% set(gca,'YTickMode','manual','YTickLabelMode','manual');

% --------------------------------------------------------------------
function exdata_Callback(hObject, eventdata, handles)


function figure1_KeyPressFcn(hObject, eventdata, handles)
aa=get(gcf,'CurrentCharacter');
set(gcf,'Pointer','arrow');
switch aa,
    case 'O', tt2dtomo('file_Callback',gcbo,[],guidata(gcbo));
    case '1', tt2dtomo('start_Callback',gcbo,[],guidata(gcbo));
    case '1', tt2dtomo('run_Callback',gcbo,[],guidata(gcbo));
    case 'C', tt2dtomo('clustermodel_Callback',gcbo,[],guidata(gcbo));
    case 'E', tt2dtomo('export_Callback',gcbo,[],guidata(gcbo));
    case 'M', tt2dtomo('showlog_Callback',gcbo,[],guidata(gcbo));
    case 'R', tt2dtomo('showray_Callback',gcbo,[],guidata(gcbo));
    case 'X', tt2dtomo('exit_Callback',gcbo,[],guidata(gcbo));
end

% --------------------------------------------------------------------
function filemenu_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function helpmenu_Callback(hObject, eventdata, handles)
set(handles.help,'Enable','Off');
set(handles.doku,'Enable','Off');

% --------------------------------------------------------------------
function help_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
msgbox({'TT2dTomo - 2d travel time tomography','version 0.3','Authors: T. Günther & C. Rücker','www.resistivity.net'},'About');

% --------------------------------------------------------------------
function doku_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function showmenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function showlog_Callback(hObject, eventdata, handles)
global Mesh velocity
mal=struct('canot','m/s','clog',1);
axes(handles.malfeld);
tripatchmod(Mesh,velocity,1,mal);

% --------------------------------------------------------------------
function showvel_Callback(hObject, eventdata, handles)
global Mesh velocity
mal=struct('canot','m/s','clog',0);
axes(handles.malfeld);
tripatchmod(Mesh,velocity,1,mal);

% --------------------------------------------------------------------
function showslow_Callback(hObject, eventdata, handles)
global Mesh velocity
mal=struct('canot','ms/m','clog',0);
axes(handles.malfeld);
tripatchmod(Mesh,1000./velocity,1,mal);

% --------------------------------------------------------------------
function showdata_Callback(hObject, eventdata, handles)
global Shot t
axes(handles.malfeld);
shotimage(Shot);

% --------------------------------------------------------------------
function showerrors_Callback(hObject, eventdata, handles)
global Shot
axes(handles.malfeld);
shotimage(Shot,0,Shot.sd./Shot.t*100);

% --------------------------------------------------------------------
function showrec_Callback(hObject, eventdata, handles)
global Shot
axes(handles.malfeld);
shotimage(Shot,1);

% --------------------------------------------------------------------
function showcov_Callback(hObject, eventdata, handles)
global Mesh W
axes(handles.malfeld);
mal=struct('canot','cov/m','clog',1);
cov=sum(W);mi=min(cov(find(cov)));
tripatchmod(Mesh,cov+mi,(cov>0),mal);


% --------------------------------------------------------------------
function showray_Callback(hObject, eventdata, handles)
global Mesh Shot velocity
axes(handles.malfeld);
showrays(Mesh,Shot,ceil(rand(1)*length(Shot.ns)),velocity);

% --------------------------------------------------------------------
function showallrays_Callback(hObject, eventdata, handles)
global Mesh Shot velocity
axes(handles.malfeld);
for i=1:length(Shot.ns), showrays(Mesh,Shot,i,velocity); end

% --------------------------------------------------------------------
function clustermodel_Callback(hObject, eventdata, handles)
global Mesh Shot t velocity W C
set(gcf,'Pointer','Watch');
velocity=clustermodel(Mesh,velocity);
mal=struct('canot','m/s');
axes(handles.malfeld);
tripatchmod(Mesh,velocity,1,mal);
set(gcf,'Pointer','Watch');drawnow;
[W,t]=waymatrix(Mesh,Shot,velocity);
chiq=chi2(Shot.t,t,0.001,0); % err
% set(handles.message,'String',...
%     sprintf('CHI^2=%.2f SQD=%.2fµs\n',chiq(end),sqrt(mean((Shot.t-t).^2))*1e6));
set(handles.message,'String',...
    sprintf('CHI^2=%.2f SQD=%ss\n',chiq(end),int2hum(sqrt(mean((Shot.t-t).^2)))));
set(gca,'xlim',[min(Mesh.node(:,1)) max(Mesh.node(:,1))]);
% coverage=sum(W);scoverage=(C'*(C*coverage(:))~=0);
scoverage=1;
axes(handles.malfeld);
tripatchmod(Mesh,velocity,scoverage,mal);
axes(handles.malfeld);
nsm=cell2mat(Shot.ns);hold on;plot(Shot.pos(nsm,1),Shot.pos(nsm,2),'k.');hold off
set(gcf,'Pointer','Arrow');
