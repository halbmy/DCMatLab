function varargout = ra2dtomo(varargin)
% RA2DTOMO M-file for ra2dtomo.fig
%      RA2DTOMO, by itself, creates a new RA2DTOMO or raises the existing
%      singleton*.
%
%      H = RA2DTOMO returns the handle to a new RA2DTOMO or the handle to
%      the existing singleton*.
%
%      RA2DTOMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RA2DTOMO.M with the given input arguments.
%
%      RA2DTOMO('Property','Value',...) creates a new RA2DTOMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ra2dtomo_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ra2dtomo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ra2dtomo

% Last Modified by GUIDE v2.5 27-Feb-2008 09:17:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ra2dtomo_OpeningFcn, ...
                   'gui_OutputFcn',  @ra2dtomo_OutputFcn, ...
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


% --- Executes just before ra2dtomo is made visible.
function ra2dtomo_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for ra2dtomo
handles.output = hObject;
set(handles.upper,'Enable','On');
% Update handles structure
guidata(hObject, handles);
set(handles.text12,'String','zweight');
set(handles.zpower,'String','0.1');
% set(handles.upper,'Enable','on');
% set(handles.dx,'String',0.5);
% UIWAIT makes ra2dtomo wait for user response (see UIRESUME)
% uiwait(handles.figure1);
% iconify(hObject,'ra2dtomo.ico');
set(gcf,'Renderer','OpenGL')

% --- Outputs from this function are returned to the command line.
function varargout = ra2dtomo_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in file.
function file_Callback(hObject, eventdata, handles)
global Shot datfile
allfiles='*.txt;*.grm;*.tom;*.dat;*.sgt;*.fbp;*.csv';
if isempty(datfile), datfile='*.dat'; end
[fname,pname]=uigetfile(allfiles,'Data File',fullfile(fileparts(datfile),allfiles));
if isequal(fname,0), return; end
datfile=fullfile(pname,fname);
[pp,ff,ee]=fileparts(datfile);
if isequal(lower(ee),'.fbp'),
    Shot=readfbpfile(datfile);    
elseif isequal(lower(ee),'.csv'),
    Shot=readcsvrafile(datfile);
elseif isequal(lower(ee),'.txt'),
    Shot=readgli(datfile);
    Shot.t=Shot.t/1000;
elseif isequal(lower(ee),'.grm'),
    Shot=readgrm(datfile);
    Shot.t=Shot.t/1000;
elseif isequal(lower(ee),'.tom'),
    Shot=readtom(datfile);
    Shot.t=Shot.t/1000;
elseif isequal(lower(ee),'.dat')|isequal(lower(ee),'.sgt'),
    fid=fopen(datfile);ss=sscanf(fgetl(fid),'%d');fclose(fid);
    if length(ss)==1, 
        Shot=readunishot(datfile);
        if size(Shot.pos,2)>2, % 3 dimensions given
%            xmbm=[0;cumsum(sqrt(sum(diff(Shot.pos(:,1:2)).^2,2)))]+Shot.pos(1,1);
%            Shot.pos(:,1)=xmbm;
           Shot.pos(:,2)=[];
        end
        if find(diff(Shot.pos)<0),
           [po,idx]=sort(Shot.pos(:,1));
           Shot.pos=Shot.pos(idx,:);
           aidx=zeros(size(idx));aidx(idx)=1:length(idx);
           for i=1:length(Shot.ns),
              Shot.ns{i}=aidx(Shot.ns{i}); 
              Shot.nx{i}=aidx(Shot.nx{i});           
           end
        end
    else Shot=readmax(datfile); end
else
    errordlg('File type not recognized!');return;
end
if ~isfield(Shot,'ns'),
    uns=unique(Shot.s);
    for i=1:length(uns),
        Shot.ns{i}=uns(i);
        fi=find(Shot.s==uns(i));
        Shot.nx{i}=Shot.g(fi);
        Shot.tt{i}=Shot.t(fi)*1000;
    end
end
Shot.pos=round(Shot.pos*1000)/1000;
axes(handles.modelaxes);cla reset;
axes(handles.dataaxes);
plotshot(Shot);
nn=[ 'Ra2DTomo - ' strrep(datfile,pwd,'')];
try, set(hObject,'Name',nn);
catch, set(handles.figure1,'Name',nn); end
ma=0;
for i=1:length(Shot.ns),
   ma1=max(abs(Shot.pos(Shot.nx{i},1)-Shot.pos(Shot.ns{i},1)));
   if ma1>ma, ma=ma1; end
end
set(handles.depth,'String',num2str(ceil(ma/3)));
set(handles.makemesh,'Enable','On');
set(handles.start,'Enable','Off');
set(handles.run,'Enable','Off');

function depth_Callback(hObject, eventdata, handles)
if ~isnumeric(str2num(get(hObject,'String'))), set(hObject,'String','0'); end


% --- Executes during object creation, after setting all properties.
function depth_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in makemesh.
function makemesh_Callback(hObject, eventdata, handles)
global Shot Mesh datfile
[pp,ff,ee]=fileparts(datfile);
meshname=strrep(datfile,ee,'');
depth=str2num(get(handles.depth,'String'));
if ~isnumeric(depth), depth=10; end
dx=str2num(get(handles.dx,'String'));
if ~isnumeric(dx), dx=0.5; end
axes(handles.dataaxes);plotshot(Shot);
axes(handles.modelaxes);cla reset;drawnow;
if 1,
%     writepoly(Shot.pos,[meshname '.poly'],depth,dx);
    Poly=pos2poly(Shot.pos,depth,dx);
    writepoly2d([meshname '.poly'],Poly);
    dos(['dctriangle -v -S -q34.3 "' meshname '.poly"']);
    Mesh=loadmesh([meshname '.bms']);
else
    Poly=create2dpoly(Shot.pos,depth,dx);
    Poly.node(:,4)=Poly.node(:,3);
    Poly.node(:,3:4)=1;
    Poly.node=unique(Shot.pos,'rows');
    Poly.node(end+1,:)=[Poly.node(end,1) Poly.node(end,2)-depth];
    Poly.node(end+1,:)=[Poly.node(1,1) Poly.node(1,2)-depth];
    Poly.node(:,3)=0.5; %triangle side
    Poly.node(:,4)=1;
    Poly.edge=(1:size(Poly.node,1))';
    Poly.edge(:,2)=Poly.edge(:,1)+1;    
    Poly.edge(end,2)=1;
    Poly.edge(:,3)=2;
    writeeapoly(Poly,[meshname '.d']);
    dos(['easymesh ' meshname '.d']);
    Mesh=readeasymesh(meshname);
%     delete([meshname '.s']);
%     delete([meshname '.e']);
%     delete([meshname '.n']);
    savemesh(Mesh,[meshname '.bms']);
end
axes(handles.modelaxes);
patch('Vertices',Mesh.node,'Faces',Mesh.cell,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]);
axis equal tight
% tripatchmod(Mesh);
% set(gca,'DataAspectRatioMode','auto','PlotBoxAspectRatioMode','auto');%daspect([1 0.5 1]);
nsm=cell2mat(Shot.ns);hold on;plot(Shot.pos(nsm,1),Shot.pos(nsm,2),'k.');hold off
set(handles.start,'Enable','On');
set(handles.run,'Enable','Off');
set(gcf,'Pointer','Arrow');

function vmin_Callback(hObject, eventdata, handles)
if ~isnumeric(str2num(get(hObject,'String'))), set(hObject,'String','0'); end


% --- Executes during object creation, after setting all properties.
function vmin_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function vmax_Callback(hObject, eventdata, handles)
if ~isnumeric(str2num(get(hObject,'String'))), set(hObject,'String','0'); end

% --- Executes during object creation, after setting all properties.
function vmax_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in start.
function start_Callback(hObject, eventdata, handles)
global Shot Mesh velocity t W chiq C offset datfile coverage
offset=[];
zweight = str2double( get(handles.zpower,'String') );
if ~isnumeric(zweight), 
    zweight = 1; 
end
C = onesidesmoothness(Mesh,zweight,1);
vmin = str2double(get(handles.vmin,'String'));
vmax = str2double(get(handles.vmax,'String'));
% velocity(:)=vmin;
% cellmids=zeros(Mesh.ncells,2);
% for i=1:Mesh.ncells, cellmids(i,:)=mean(Mesh.node(Mesh.cell(i,:),:)); end
% mi=min(Mesh.node(:,2));ma=max(Mesh.node(:,2));
% velocity(cellmids(:,2)<mi+(ma-mi)/3)=vmax;
velocity=getstartmodel(Mesh,vmin,vmax);
mal=struct('canot','m/s');
axes(handles.modelaxes);tripatchmod(Mesh,velocity,1,mal);
set(gcf,'Pointer','Watch');
axes(handles.dataaxes);cla reset;drawnow;
[W,t]=waymatrix(Mesh,Shot,velocity);
chiq=chi2(Shot.t,t,0.001,0); % err
out=sprintf(' CHI^2=%.2f RMS=%.2fms\n',chiq(end),sqrt(mean((Shot.t-t).^2))*1000);
nn=[ 'Ra2DTomo - ' strrep(datfile,pwd,'') ' -' out];
try, set(hObject,'Name',nn);
catch, set(handles.figure1,'Name',nn); end
plotshot(Shot,t*1.000);
set(gca,'xlim',[min(Mesh.node(:,1)) max(Mesh.node(:,1))]);
coverage=sum(W);
scoverage=(C'*(C*coverage(:))~=0);
axes(handles.modelaxes);tripatchmod(Mesh,velocity,scoverage,mal);
% set(gca,'DataAspectRatioMode','auto','PlotBoxAspectRatioMode','auto');%daspect([1 0.5 1]);
nsm=cell2mat(Shot.ns);hold on;plot(Shot.pos(nsm,1),Shot.pos(nsm,2),'k.');hold off
axes(handles.modelaxes);
%showrays(Mesh,Shot,fix(length(Shot.ns)/2),velocity);
set(handles.run,'Enable','On');
set(gcf,'Pointer','Arrow');


function lower_Callback(hObject, eventdata, handles)
if ~isnumeric(str2num(get(hObject,'String'))), set(hObject,'String','0'); end

% --- Executes during object creation, after setting all properties.
function lower_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function upper_Callback(hObject, eventdata, handles)
if ~isnumeric(str2num(get(hObject,'String'))), set(hObject,'String','0'); end

% --- Executes during object creation, after setting all properties.
function upper_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in parameter.
function parameter_Callback(hObject, eventdata, handles)
pa=get(handles.parameter,'Value');
onoff='Off';
if pa==1, onoff='On'; end
set(handles.lower,'Enable',onoff);
set(handles.upper,'Enable',onoff);

% --- Executes during object creation, after setting all properties.
function parameter_CreateFcn(hObject, eventdata, handles)
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
global Mesh velocity Shot t W C chiq offset datfile coverage
try,
withoffset=1;
nModel=length(velocity);nShot=length(Shot.ns);
if withoffset,
    if isempty(offset)||(length(offset)~=nShot), offset=zeros(nShot,1); end
    if size(W,2)==nModel, [W,t]=waymatrix(Mesh,Shot,velocity,[],offset); end
else
    offset=[];
end
deltaData=Shot.t-t;err=1/1000;
D=1/err;D2=D*D;lv=length(velocity);
lbound=str2num(get(handles.lower,'String'));
if lbound>min(velocity), lbound=0;set(handles.lower,'String','0'); end
ubound=str2num(get(handles.upper,'String'));
if ubound<max(velocity), ubound=0;set(handles.upper,'String','0'); end
lam=str2double(get(handles.lambda,'String'));
par=get(handles.parameter,'Value');
axes(handles.modelaxes);cla reset;drawnow;
set(gcf,'Pointer','Watch');
switch par, % Inversion parameter
    case 1, % log(velocity)
        deltaModel=logtrans(velocity,lbound,ubound);
        vv=-1./(velocity.^2.*logtransdiff(velocity,lbound,ubound));
    case 2, % velocity
        deltaModel=velocity;
        vv=-1./velocity.^2;lam=lam/1e6;
    case 3, % slowness
        deltaModel=1./velocity;
        vv=ones(lv,1);lam=lam*1e6;
end
if withoffset, 
    vv=[vv;ones(nShot,1)]; 
%     deltaModel=[deltaModel;zeros(length(offset),1)]; %local
    deltaModel=[deltaModel;offset]; %global
    olam=lam*2;
    nBound=size(C,1);
    for i=1:nShot, C(nBound+i,nModel+i)=olam; end
end
S=W*spdiags(vv,0,length(vv),length(vv));
% ss=svd(full(S));max(ss)
%     dm=cglscdp(S,deltaData,lam,C,D,1,deltaModel);
dm=(S'*D2*S+lam*C'*C)\(S'*(D2*deltaData)-lam*(C'*C*deltaModel));
if withoffset,
    offset=offset+dm(nModel+1:end);
    fprintf('Offsets:');fprintf(' %.1f',offset*1000);fprintf('ms\n');
    dm(nModel+1:end)=[];
    C(nBound+1:end,:)=[];
    C(:,nModel+1:end)=[];
end
switch par,
    case 1, velocity1=ilogtrans(logtrans(velocity,lbound,ubound)+dm,lbound,ubound);
    case 2, velocity1=velocity+dm;
    case 3, velocity1=1./(1./velocity+dm);
end
fprintf('min/max v=%d/%d',round(min(velocity1)),round(max(velocity1)));
mal=struct('canot','m/s','clog',1,'cauto',1);
tripatchmod(Mesh,velocity1,1,mal);
axes(handles.dataaxes);cla reset;drawnow;
oldt=t;[W,t]=waymatrix(Mesh,Shot,velocity1,[],offset);
for i=1:21, % line search procedure
    tau=0.05*i;ta=tau*t+(1-tau)*oldt;
    phia(i)=chi2(Shot.t,ta,err,0);
end
[pp,ii]=min(phia);fak=ii*0.05;
if fak<1,
    fprintf('  Line search factor %.2f',fak);dm=dm*fak;
    switch par,
        case 1,
            if ubound>0,
                velocity=(ubound*(velocity-lbound).*exp(dm)+lbound*(ubound-velocity))./((velocity-lbound).*exp(dm)+ubound-velocity);
            else
                velocity=(velocity-lbound).*exp(dm)+lbound;
            end
        case 2, velocity=velocity+dm;
        case 3, velocity=1./(1./velocity+dm);
    end
    fprintf('  min/max v=%d/%d',round(min(velocity)),round(max(velocity)));
    [W,t]=waymatrix(Mesh,Shot,velocity,[],offset);
else
    velocity=velocity1;
end
coverage=sum(W);
scoverage=(C'*(C*coverage(1:nModel)')~=0);
axes(handles.modelaxes);tripatchmod(Mesh,velocity,scoverage,mal);%axis normal
nsm=cell2mat(Shot.ns);hold on;
plot(Shot.pos(nsm,1),Shot.pos(nsm,2),'k.');hold off
% set(gca,'DataAspectRatioMode','auto','PlotBoxAspectRatioMode','auto');%daspect([1 0.5 1]);
axes(handles.dataaxes);plotshot(Shot,t*1.000,offset);
set(gca,'xlim',[min(Mesh.node(:,1)) max(Mesh.node(:,1))]);
chiq(end+1)=chi2(Shot.t,t,err,0);
out=sprintf(' CHI^2=%.2f RMS=%.2fms\n',chiq(end),sqrt(mean((Shot.t-t).^2))*1000);
fprintf(out);
nn=[ 'Ra2DTomo - ' strrep(datfile,pwd,'') ' -' out];
try, set(hObject,'Name',nn);
catch, set(handles.figure1,'Name',nn); end
if get(handles.blocky,'Value'), % blocky inversion
    dx=str2double(get(handles.dx,'String'));
    C1=onesidesmoothness(Mesh,dx);
    sm=abs(C1*log(velocity));su2=sum(sm.^2);sua=sum(sm);
    wxz=ones(size(sm));fi=find(sm);
    if ~isempty(fi), wxz(fi)=su2/sua./sm(fi); end
    wxz(wxz>1)=1;
    C=spdiags(wxz,0,length(wxz),length(wxz))*C1;
end
axes(handles.modelaxes);
% showrays(Mesh,Shot,ceil(rand(1)*length(Shot.ns)),velocity);
set(gcf,'Pointer','Arrow');
catch,
    set(gcf,'Pointer','Arrow');
    disp(lasterr);
end

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
if ~isnumeric(str2double(get(hObject,'String'))), 
    set(hObject,'String','0'); 
end
zweight = str2double(get(handles.zpower,'String'));
C=onesidesmoothness(Mesh,zweight,1);

% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)
global datfile
[pp,nn,ee]=fileparts(datfile);
[fname,pname]=uiputfile('*.eps','Export Figure',strrep(datfile,ee,'_fig.eps'));
if ischar(fname),
%     print(handles.figure1,'-depsc2',fullfile(pname,fname));
    epsprint(handles.figure1,fullfile(pname,fname));
end

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
ra2dtomo('file_Callback',gcbo,[],guidata(gcbo))

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
global Mesh velocity Shot
onoff='off';
if isfield(Mesh,'ncells')&&(Mesh.ncells==length(velocity)), onoff='On'; end
set(handles.exmodel,'Enable',onoff);
onoff='off';
if isstruct(Shot), onoff='on'; end
set(handles.exdata,'Enable',onoff);

% --------------------------------------------------------------------
function exfigure_Callback(hObject, eventdata, handles)
ra2dtomo('export_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function exmodel_Callback(hObject, eventdata, handles)
global Mesh velocity W C datfile coverage
f=figure(1);
set(1,'MenuBar','none','NumberTitle','off','Name','Data,Model response & misfit');
% iconify(1,'ra2dromo.ico'); 
% coverage=sum(W)';
scoverage=(C'*(C*coverage(1:Mesh.ncells))~=0);
mal=struct('canot','m/s','clog',1,'cauto',1);
set(handles.modelaxes,'Units','pixel');
po=get(handles.modelaxes,'Position');
po(1:2)=[20 40];po(4)=po(4)+100;
set(f,'Units','pixel','position',po); % gleiche figure-groesse
set(handles.modelaxes,'Units','normalized'); % zurück
set(f,'PaperSize',po(3:4),'PaperPositionMode','auto');
tripatchmod(Mesh,velocity,scoverage,mal);
% set(gca,'YTickMode','manual','YTickLabelMode','manual');
[pp,nn,ee]=fileparts(datfile);
[fname,pname]=uiputfile('*.eps','Export Graphics',strrep(datfile,ee,'.eps'));
if ischar(fname),
%     print(handles.figure1,'-depsc2',fullfile(pname,fname));
    epsprint(f,fullfile(pname,fname));
end
close(f);

% --------------------------------------------------------------------
function exdata_Callback(hObject, eventdata, handles)
global Shot datfile
f=figure(1);
set(1,'MenuBar','none','NumberTitle','off','Name','Data,Model response & misfit');
% iconify(1,'ra2dromo.ico'); 
set(handles.modelaxes,'Units','pixel');
po=get(handles.modelaxes,'Position');
po(1:2)=[20 40];po(4)=po(4)+150;
set(f,'Units','pixel','position',po); % gleiche figure-groesse
set(handles.dataaxes,'Units','normalized'); % zurück
set(f,'PaperSize',po(3:4),'PaperPositionMode','auto');
plotshot(Shot)
% set(gca,'YTickMode','manual','YTickLabelMode','manual');
[pp,nn,ee]=fileparts(datfile);
[fname,pname]=uiputfile('*.eps','Export Graphics',strrep(datfile,ee,'_data.eps'));
if ischar(fname),
    epsprint(f,fullfile(pname,fname));
end
close(f);


function figure1_KeyPressFcn(hObject, eventdata, handles)
aa=get(gcf,'CurrentCharacter');
set(gcf,'Pointer','arrow');
switch aa,
    case 'O', ra2dtomo('file_Callback',gcbo,[],guidata(gcbo));
    case '0', ra2dtomo('start_Callback',gcbo,[],guidata(gcbo));
    case '1', ra2dtomo('run_Callback',gcbo,[],guidata(gcbo));
    case 'C', ra2dtomo('clustermodel_Callback',gcbo,[],guidata(gcbo));
    case 'E', ra2dtomo('export_Callback',gcbo,[],guidata(gcbo));
    case 'L', ra2dtomo('showlog_Callback',gcbo,[],guidata(gcbo));
    case 'M', ra2dtomo('makemesh_Callback',gcbo,[],guidata(gcbo));
    case 'R', ra2dtomo('showray_Callback',gcbo,[],guidata(gcbo));
    case 'S', ra2dtomo('showslow_Callback',gcbo,[],guidata(gcbo));
    case 'V', ra2dtomo('showvel_Callback',gcbo,[],guidata(gcbo));
    case 'X', ra2dtomo('exit_Callback',gcbo,[],guidata(gcbo));
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
msgbox({'Ra2dTomo - 2d Refraction tomography','version 0.8.2','Authors: T. Günther & C. Rücker','www.resistivity.net'},'About');

% --------------------------------------------------------------------
function doku_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function showmenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function showlog_Callback(hObject, eventdata, handles)
global Mesh velocity W C coverage
% coverage=sum(W)';
scoverage=(C'*(C*coverage(1:Mesh.ncells))~=0);
mal=struct('canot','m/s','clog',1,'cauto',1);
axes(handles.modelaxes);
tripatchmod(Mesh,velocity,scoverage,mal);

% --------------------------------------------------------------------
function showvel_Callback(hObject, eventdata, handles)
global Mesh velocity W C coverage
% coverage=sum(W)';
if length(coverage)==Mesh.ncells,
    scoverage=(C'*(C*coverage(1:Mesh.ncells))~=0);
else
    scoverage=ones(Mesh.ncells,1);
end
mal=struct('canot','m/s','clog',0,'cauto',1);
axes(handles.modelaxes);tripatchmod(Mesh,velocity,scoverage,mal);

% --------------------------------------------------------------------
function showslow_Callback(hObject, eventdata, handles)
global Mesh velocity W C coverage
% coverage=sum(W)';
if length(coverage)==Mesh.ncells,
    scoverage=(C'*(C*coverage(1:Mesh.ncells))~=0);
else
    scoverage=ones(Mesh.ncells,1);
end
mal=struct('canot','ms/m','clog',0,'cauto',1);
axes(handles.modelaxes);tripatchmod(Mesh,1000./velocity,scoverage,mal);

% --------------------------------------------------------------------
function showdata_Callback(hObject, eventdata, handles)
global Shot t
axes(handles.dataaxes);
plotshot(Shot,t*1.000);

% --------------------------------------------------------------------
function showshotimage_Callback(hObject, eventdata, handles)
global Shot
figure(1);
set(1,'Units','Characters','MenuBar','none','NumberTitle','off','Name','Shot image');
% iconify(1,'ra2dtomo.ico');
% shotimage(Shot);
va=getva(Shot);
showasimage(Shot,va,[],Shot.pos(:,1));
% showasimage(Shot,Shot.t);

% --------------------------------------------------------------------
function showrec_Callback(hObject, eventdata, handles)
global Shot
figure(1);
set(1,'Units','Characters','MenuBar','none','NumberTitle','off','Name','Shot image');
% iconify(1,'ra2dtomo.ico');
shotimage(Shot,1);

% --------------------------------------------------------------------
function showcov_Callback(hObject, eventdata, handles)
global Mesh W coverage
axes(handles.modelaxes);
mal=struct('canot','cov/m','clog',1,'cauto',1);
% cov=sum(W);
mi=min(coverage(find(coverage)));
tripatchmod(Mesh,coverage+mi,(coverage>0),mal);


% --------------------------------------------------------------------
function showray_Callback(hObject, eventdata, handles)
global Mesh Shot velocity
axes(handles.modelaxes);
showrays(Mesh,Shot,ceil(rand(1)*length(Shot.ns)),velocity);

% --------------------------------------------------------------------
function showallrays_Callback(hObject, eventdata, handles)
global Mesh Shot velocity
axes(handles.modelaxes);
for i=1:length(Shot.ns), showrays(Mesh,Shot,i,velocity); end

% --------------------------------------------------------------------
function clustermodel_Callback(hObject, eventdata, handles)
global Mesh Shot t velocity W C coverage
set(gcf,'Pointer','Watch');
velocity=clustermodel(Mesh,velocity);
mal=struct('canot','m/s');
axes(handles.modelaxes);tripatchmod(Mesh,velocity,1,mal);
set(gcf,'Pointer','Watch');drawnow;
[W,t]=waymatrix(Mesh,Shot,velocity);
chiq=chi2(Shot.t,t,0.001,0); % err
fprintf('CHI^2=%.2f RMS=%.2fms\n',chiq(end),sqrt(mean((Shot.t-t).^2))*1000);
axes(handles.dataaxes);plotshot(Shot,t*1.000);
set(gca,'xlim',[min(Mesh.node(:,1)) max(Mesh.node(:,1))]);
nModel=length(velocity);
% coverage=sum(W);
scoverage=(C'*(C*coverage(1:nModel)')~=0);
axes(handles.modelaxes);tripatchmod(Mesh,velocity,scoverage,mal);
nsm=cell2mat(Shot.ns);hold on;plot(Shot.pos(nsm,1),Shot.pos(nsm,2),'k.');hold off
% axes(handles.modelaxes);showrays(Mesh,Shot,fix(length(Shot.ns)/2),velocity);
set(gcf,'Pointer','Arrow');


% --------------------------------------------------------------------
function runmenu_Callback(hObject, eventdata, handles)
% hObject    handle to runmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function runstep_Callback(hObject, eventdata, handles)
ra2dtomo('run_Callback',gcbo,[],guidata(gcbo))


% --------------------------------------------------------------------
function runfull_Callback(hObject, eventdata, handles)
global Mesh velocity Shot t W C chiq offset
running=1;
while running,
    ra2dtomo('run_Callback',gcbo,[],guidata(gcbo))
    running=(chiq(end)/chiq(end-1)<0.95)&(chiq(end)>1);
end
    

% --------------------------------------------------------------------
function runttinv_Callback(hObject, eventdata, handles)
global Mesh velocity Shot t W C chiq offset datfile coverage
[pp,ff,ee]=fileparts(datfile);
meshname=strrep(datfile,ee,'.bms');
if ismember(ee,{'.txt','.grm','.tom'})
    datfile=fullfile(pp,[ff '.sgt']);
    savesgtfile(datfile,Shot);
end
lbound = str2double(get(handles.lower,'String'));
ubound = str2double(get(handles.upper,'String'));
lam = str2double(get(handles.lambda,'String'));
zweight = str2double(get(handles.zpower,'String'));
err = 1/1000;
scall=['ttinv -vG -p ' meshname ' -e 0.1 -l ' num2str(lam) ' -t ' num2str(err) ' '];
if lbound>0, scall=[scall '-u ' num2str(1/lbound) ' ']; end
if ubound>lbound, scall=[scall '-b ' num2str(1/ubound) ' ']; end
if zpower>0, scall=[scall '-z ' num2str(zweight) ' ']; end
if get(handles.blocky,'Value'), scall=[scall '-B ']; end
if get(handles.robust,'Value'), scall=[scall '-R ']; end
scall = [scall datfile];
fprintf('%s\n',scall);
systemcall(scall);
fid=fopen('velocity.vector');velocity=fscanf(fid,'%f',Mesh.ncells);fclose(fid);
fid=fopen('coverage.vector');coverage=fscanf(fid,'%f',Mesh.ncells);fclose(fid);
fid=fopen('response.vector');t=fscanf(fid,'%f',length(Shot.t));fclose(fid);
mal=struct('canot','m/s');
if (nnz(C)==0)|(Mesh.ncells~=size(C,2)), C=onesidesmoothness(Mesh); end
scoverage=(C'*(C*coverage(:))~=0);
axes(handles.modelaxes);tripatchmod(Mesh,velocity,scoverage,mal);
axes(handles.dataaxes);plotshot(Shot,t);%*1000);
scoverage=(scoverage+0.1)/1.1;
save('scoverage.vector','scoverage','-ascii');

% --------------------------------------------------------------------
function modelascii_Callback(hObject, eventdata, handles)
global Mesh velocity datfile
if Mesh.ncells~=length(velocity), error('Model size mismatch!'); end
xzv=meshcellmid(Mesh);xzv(:,2)=-xzv(:,2);
xzv(:,3)=velocity(:);
[pp,nn,ee]=fileparts(datfile);
[fname,pname]=uiputfile('*.xzv','Export Model as Ascii',strrep(datfile,ee,'.xzv'));
if ischar(fname),
    fid=fopen(fullfile(pname,fname),'w');
    fprintf(fid,'x[m]\tz[m]\tv[m/s]\r\n');
    fprintf(fid,'%.2f\t%.2f\t%.2f\r\n',xzv');
    fclose(fid);
end