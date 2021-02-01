function varargout = underwater(varargin)
% UNDERWATER M-file for underwater.fig
%      UNDERWATER, by itself, creates a new UNDERWATER or raises the existing
%      singleton*.
%
%      H = UNDERWATER returns the handle to a new UNDERWATER or the handle to
%      the existing singleton*.
%
%      UNDERWATER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNDERWATER.M with the given input arguments.
%
%      UNDERWATER('Property','Value',...) creates a new UNDERWATER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before underwater_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to underwater_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help underwater

% Last Modified by GUIDE v2.5 22-Feb-2008 14:39:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @underwater_OpeningFcn, ...
                   'gui_OutputFcn',  @underwater_OutputFcn, ...
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


% --- Executes just before underwater is made visible.
function underwater_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to underwater (see VARARGIN)

% Choose default command line output for underwater
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes underwater wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = underwater_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editdatafile_Callback(hObject, eventdata, handles)
global DATA datafile
datafile=get(handles.editdatafile,'String');
ilon=0;ilat=0;ii=0;iu=0;ip=0;id=0;itime=0;
%%
fid=fopen(datafile);
if fid<0, 
    set(handles.statusdatafile,'String',['Could not open ' datafile]);
    return;
end
zeile=fgetl(fid);%fclose(fid);
i=0;formstr='';
while ~isempty(zeile),
    i=i+1;
    [tok,zeile]=strtok(zeile);ss='%f';
    if strcmp(lower(tok),'lon'), ilon=i;ss='%f %*s'; 
    elseif strcmp(lower(tok),'lat'), ilat=i;ss='%f %*s';
    elseif strcmp(lower(tok),'i'), ii=i;
    elseif strcmp(lower(tok),'u'), iu=i;
    elseif strcmp(lower(tok),'p'), ip=i;
    elseif strcmp(lower(tok),'d'), id=i;
    elseif strcmp(lower(tok),'channel'), ich=i;ss='%d';
    elseif strcmp(lower(tok),'time'), itime=i;ss='%s %*s';
    else ss='%*s';i=i-1; end
    formstr=[formstr ss ' '];
end
formstr(end)='';
[All{1},All{2},All{3},All{4},All{5},All{6},All{7},All{8}]=textread(datafile,formstr,'headerlines',1);
% All=textscan(fid,formstr);
fclose(fid);
DATA=[];DATA.ndata=0;
if iu, DATA.u=All{iu};DATA.ndata=length(DATA.u); end
if ii, DATA.i=All{ii}; end
if id, DATA.d=All{id}; end
if ip, DATA.p=All{ip}; end
if ilat, DATA.lat=All{ilat}; end
if ilon, DATA.lon=All{ilon}; end
if itime, 
    DATA.tstr=All{itime}; 
    t0=datenum('00:00:00');
    DATA.time=datenum(DATA.tstr)-t0;
end
set(handles.statusdatafile,'String',sprintf('%d data (%s-%s)',...
    DATA.ndata,datestr(min(DATA.time),15),datestr(max(DATA.time),15)));
%%

% --- Executes during object creation, after setting all properties.
function editdatafile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in browsedatafile.
function browsedatafile_Callback(hObject, eventdata, handles)
global datafile
filetype='*.txt';infile=filetype;
if ~isempty(datafile),
    [pp,ff,ee]=fileparts(datafile);
    infile=fullfile(pp,infile);
end
[fname,pname]=uigetfile(filetype,'Load data file',infile);
if ~isstr(fname), return; end
datafile=fullfile(pname,fname);
set(handles.editdatafile,'String',datafile);
underwater('editdatafile_Callback',gcbo,[],guidata(gcbo));


function editgeolorefile_Callback(hObject, eventdata, handles)
global Geolore
geolorefile=get(handles.editgeolorefile,'String');
A=load(geolorefile);
Geolore=[];
Geolore.time=A(:,2)/24;
Geolore.depth=(A(:,3)-0.456)/0.043;
gpsfile=get(handles.editgpsfile,'String');
if exist(gpsfile,'file')
    dd=textread(gpsfile,'%*s%f%*s',1,'headerlines',3);
    nulltime=floor(dd/10000)/24+floor(mod(dd,10000)/100)/24/60+mod(dd,100)/86400;
    Geolore.time=Geolore.time+nulltime;
end
set(handles.statusgeolorefile,'String',sprintf('%d readings (%s-%s), d=%.1f-%.1f',...
    length(Geolore.time),datestr(min(Geolore.time),15),datestr(max(Geolore.time),15),min(Geolore.depth),max(Geolore.depth)));

% --- Executes during object creation, after setting all properties.
function editgeolorefile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in browsegeolorefile.
function browsegeolorefile_Callback(hObject, eventdata, handles)
global datafile
filetype='*.txt';infile=filetype;
if ~isempty(datafile),
    [pp,ff,ee]=fileparts(datafile);
    infile=fullfile(pp,infile);
end
[fname,pname]=uigetfile(filetype,'Load Geolore file',infile);
if ~isstr(fname), return; end
geolorefile=fullfile(pname,fname);
set(handles.editgeolorefile,'String',geolorefile);
underwater('editgeolorefile_Callback',gcbo,[],guidata(gcbo));



function editgpsfile_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function editgpsfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browsegpsfile.
function browsegpsfile_Callback(hObject, eventdata, handles)
global datafile
filetype='*.gps';infile=filetype;
if ~isempty(datafile),
    [pp,ff,ee]=fileparts(datafile);
    infile=fullfile(pp,infile);
end
[fname,pname]=uigetfile(filetype,'Load Geolore GPS file',infile);
if ~isstr(fname), return; end
gpsfile=fullfile(pname,fname);
set(handles.editgpsfile,'String',gpsfile);
underwater('editgeolorefile_Callback',gcbo,[],guidata(gcbo));


function editdepthlogfile_Callback(hObject, eventdata, handles)
global Log
depthlogfile=get(handles.editdepthlogfile,'String');
[HH,MM,SS,HU,D15,D100]=textread(depthlogfile,'%d:%d:%d.%d,%d,%d\r\n%*c');
Log.time=(((SS+HU/100)/60+MM)/60+HH-2)/24;
Log.depth=D15/100;
Log.depth2=D100/100;
set(handles.statusdepthlogfile,'String',sprintf('%d readings (%s-%s), d=%.1f-%.1f',...
    length(Log.time),datestr(min(Log.time),15),datestr(max(Log.time),15),min(Log.depth),max(Log.depth)));

% --- Executes during object creation, after setting all properties.
function editdepthlogfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browsedepthlogfile.
function browsedepthlogfile_Callback(hObject, eventdata, handles)
global datafile
filetype='*.txt';infile=filetype;
if ~isempty(datafile),
    [pp,ff,ee]=fileparts(datafile);
    infile=fullfile(pp,infile);
end
[fname,pname]=uigetfile(filetype,'Load depth log file',infile);
if ~isstr(fname), return; end
gpslogfile=fullfile(pname,fname);
set(handles.editdepthlogfile,'String',gpslogfile);
underwater('editdepthlogfile_Callback',gcbo,[],guidata(gcbo));



function editgpslogfile_Callback(hObject, eventdata, handles)
global Geolore
gpslogfile=get(handles.editgpslogfile,'String');
[S1,S2,datum,zeit,Geolore.lat,Geolore.lon]=textread(gpslogfile,'%d/%d%d%d%*s%f%*s%f');
Geolore.gpstime=floor(zeit/10000)/24+floor(mod(zeit,10000)/100)/24/60+mod(zeit,100)/86400;
if 1,
    [Geolore.N,Geolore.E]=gps2xyz(Geolore.lat,Geolore.lon);
    dl=sqrt(diff(Geolore.N).^2+diff(Geolore.E).^2);
    Geolore.v=dl./diff(Geolore.gpstime)/24/1000; %kn
else
    dl=sqrt(diff(Geolore.lat).^2+diff(Geolore.lon).^2);
    Geolore.v=dl./(diff(Geolore.gpstime)*24);
end
Geolore.x=cumsum([0;dl])*1000;
% minmax(Geolore.v)
set(handles.statusgpslogfile,'String',sprintf('%d readings (%s-%s) %.1fx%.1fkm',...
    length(Geolore.gpstime),datestr(min(Geolore.gpstime),15),datestr(max(Geolore.gpstime),15),...
    diff(minmax(Geolore.lat)),diff(minmax(Geolore.lon))));


% --- Executes during object creation, after setting all properties.
function editgpslogfile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browsegpslogfile.
function browsegpslogfile_Callback(hObject, eventdata, handles)
global datafile
filetype='*.txt';infile=filetype;
if ~isempty(datafile),
    [pp,ff,ee]=fileparts(datafile);
    infile=fullfile(pp,infile);
end
[fname,pname]=uigetfile(filetype,'Load gps log file',infile);
if ~isstr(fname), return; end
gpslogfile=fullfile(pname,fname);
set(handles.editgpslogfile,'String',gpslogfile);
underwater('editgpslogfile_Callback',gcbo,[],guidata(gcbo));


% --- Executes on button press in showroute.
function showroute_Callback(hObject, eventdata, handles)
global DATA Geolore
figure(11);
%%
tt=linspace(min(DATA.time),max(DATA.time),20);
if isfield(Geolore,'N'),
    [y,x]=gps2xyz(DATA.lat,DATA.lon);
    plot(x,y,'b.');
    axis equal;grid on;
    hold on;plot(Geolore.E,Geolore.N,'r.');hold off;
    tx=interp1(DATA.time(1:16:end),x(1:16:end),tt);
    ty=interp1(DATA.time(1:16:end),y(1:16:end),tt);
elseif isfield(Geolore,'lat'),
    plot(DATA.lon,DATA.lat,'b.');
    axis equal;grid on
    hold on;plot(Geolore.lon,Geolore.lat,'r.');hold off;
    tx=interp1(DATA.time(1:16:end),DATA.lon(1:16:end),tt);
    ty=interp1(DATA.time(1:16:end),DATA.lat(1:16:end),tt);
end
for i=1:length(tt),    
    set(text(tx(i),ty(i),datestr(tt(i),15)),'HorizontalAlignment','left');
end
%%

% --- Executes on button press in showdepths.
function showdepths_Callback(hObject, eventdata, handles)
figure(12);clf;
%%
global Geolore Log
if isfield(Geolore,'v'),
   subplot(2,1,1);
   plot(Geolore.gpstime(1:end-1)+diff(Geolore.gpstime)/2,Geolore.v/1.85,'b-');
   axis xy;datetick;grid on;
   xlim(minmax(Geolore.time));
   xlabel('time');ylabel('knots');
   subplot(2,1,2);
end
plot(Geolore.time,Geolore.depth,'b-');
ylim([0 max(Geolore.depth)]);
if isfield(Log,'depth'),
    hold on;plot(Log.time,Log.depth,'r-');hold off;
    ylim([0 max(Log.depth)]);
    legend('sensor','log');
else
    title('sensor depth');
end
axis ij;xlabel('time');ylabel('y in m');datetick;grid on;
xlim(minmax(Geolore.time));
%%

% --- Executes on button press in showrawdata.
function showrawdata_Callback(hObject, eventdata, handles)
global DATA Geolore
%%
DATA.nsond=fix(length(DATA.u)/16);
DATA.R=zeros(DATA.nsond,16);
DATA.Err=zeros(size(DATA.R));
for i=1:DATA.nsond,
    ii=(1:16)+(i-1)*16;
    DATA.R(i,:)=DATA.u(ii)'./DATA.i(ii)';
    DATA.Err(i,:)=DATA.d(ii)';
end
DATA.Sond=readsnd2d('sond.s2d');
DATA.Sond.elec(:,2)=5;
if isfield(Geolore,'depth'), DATA.Sond.elec(:,2)=median(Geolore.depth); end
DATA.Sond.k=getkonf2d(DATA.Sond);
DATA.Rhoa0=DATA.R*diag(DATA.Sond.k);
DATA.T=DATA.time(1:16:end);
DATA.Rhoa0(DATA.Rhoa0<0)=NaN;
DATA.Rhoa0(DATA.Rhoa0>interperc(DATA.Rhoa0,99.9))=NaN;
DATA.Rhoa0(DATA.Err>20)=NaN;
DATA.Rhoa0(:,4)=NaN;
DATA.R(isnan(DATA.Rhoa0))=NaN;
DATA.cax=10.^interperc(log10(DATA.Rhoa0),[2 98]);
DATA.cax=[1 3.5];
figure(1);imageline(DATA.T,DATA.Rhoa0,DATA.cax);



function editdof_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function editdof_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in processdata.
function processdata_Callback(hObject, eventdata, handles)
global DATA Data Geolore
nc=str2num(get(handles.editdof,'String'));
Data=[];
%%
dt=15/86400;
Data.T=(min(DATA.T):dt:max(DATA.T))';
% besser äquidistant aus den GPS
Data.R=zeros(length(Data.T),16);Data.Rhoa0=Data.R;
for i=1:16,
%     spalte=DATA.R(:,i);
    spalte=log10(DATA.Rhoa0(:,i));
    mid=spalte;mid(2:end-1)=mid(2:end-1)/2+(spalte(1:end-2)+spalte(3:end))/4;
    spalte(abs(mid./spalte-1)>0.2)=NaN;%outlier
    fi=find(isfinite(spalte));
    figure(2);
    if ~isempty(fi),
%         Data.R(:,i)=harmfit(DATA.T(fi),spalte(fi),nc,Data.T);        
        Data.Rhoa0(:,i)=10.^harmfit(DATA.T(fi),spalte(fi),nc,Data.T);
        Data.R(:,i)=Data.Rhoa0(:,i)./DATA.Sond.k(i);
%         plot(DATA.T,spalte,'b-',Data.T,Data.Rhoa0(:,i),'r-');pause;
    end
end
Data.Rhoa0(Data.Rhoa0<=0)=NaN;
if isfield(DATA,'cax'), Data.cax=DATA.cax;
else Data.cax=interperc(Data.Rhoa0,[2 98]); end
figure(2);imageline(Data.T,Data.Rhoa0,Data.cax);
dnew=harmfit(Geolore.time,Geolore.depth,30,Data.T);
Data.Rhoa=zeros(size(Data.R));
Sond=DATA.Sond;
for i=1:length(Data.T),
    Sond.elec(:,2)=dnew(i);
    knew=getkonf2d(Sond);
    Data.Rhoa(i,:)=Data.R(i,:).*knew';
end
Data.Rhoa(:,4)=NaN;
Data.Rhoa(Data.Rhoa<=0)=NaN;
figure(3);imageline(Data.T,Data.Rhoa,Data.cax);
set(handles.editmintime,'String',datestr(min(Data.T),15));
set(handles.editmaxtime,'String',datestr(max(Data.T),15));


% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)
global datafile
names={'raw','rhoa0','rhoa','','','','','','','','route','depth'};
[pname,fname,ename]=fileparts(datafile);
outname=fullfile(pname,fname);
for i=1:length(names),
    if ishandle(i), epsprint(i,[outname '-' names{i}],1); end
end



function editmintime_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function editmintime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editmaxtime_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function editmaxtime_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in writedata.
function writedata_Callback(hObject, eventdata, handles)
global Data DATA Geolore Log datafile Poly
null=datenum('00:00');
tmin=datenum(get(handles.editmintime,'String'),15)-null;
tmax=datenum(get(handles.editmaxtime,'String'),15)-null;
fi=find((Data.T>=tmin)&(Data.T<=tmax));
dx=sqrt(diff(DATA.lat(1:16:end)).^2+diff(DATA.lon(1:16:end)).^2);
% DATA.x=[0;cumsum(dx)]*1000;
% pos=interp1(DATA.T,DATA.x,Data.T(fi));
pos=interp1(Geolore.gpstime,Geolore.x,Data.T(fi));
%%
ch=1:16;ch(4)=[];
Sond=rmfield(delmeasurement(DATA.Sond,4),{'r','k','eind','names','nr'});
Sond.elec(:,2)=0;
du=Sond.a;Sond.a=Sond.m;Sond.m=du;
du=Sond.b;Sond.b=Sond.n;Sond.n=du;
oldpos=pos(1);
for i=1:length(pos),
    V=Sond;
    V.elec(:,1)=V.elec(:,1)+pos(i);
    V.rho=Data.R(i,ch);
    if i==1, 
        Out=V; 
    else
        if pos(i)-oldpos>40, Out=combdata2d(Out,V);oldpos=pos(i); end
    end
end
dep=interp1(Geolore.time,Geolore.depth,Geolore.gpstime);
Out.elec(:,2)=interp1(Geolore.x,dep,Out.elec(:,1),'linear','extrap');
[pname,fname,ename]=fileparts(datafile);
newdir=[strrep(get(handles.editmintime,'String'),':','') '-' strrep(get(handles.editmaxtime,'String'),':','')];
oldpwd=pwd;
cd(pname);
mkdir(newdir);
cd(oldpwd);
pname=[pname filesep newdir];
% mkdir(pname);
outfile=fullfile(pname,fname);
errVolt=1e-3;current=2.0;
Out.err=errVolt/current./Out.rho+0.03;
if min(diff(Out.elec(:,1)))<=0, error('Electrodes have to be sorted!'); end
saveinv2dfile([outfile '.dat'],Out);
Out.elec(:,2)=-Out.elec(:,2);
saveinv2dfile([outfile '.ohm'],Out);
% fi2=find((Log.time>=tmin)&(Log.time<=tmax));
de=interp1(Log.time+mod(datenum('0:0:25'),1),Log.depth,Data.T(fi),'linear','extrap');
plot(pos,de,'r-',Out.elec(:,1),Out.elec(:,2),'b-');axis ij;legend('log','sensor');
xz=[pos de];save([outfile '.xd'],'xz','-ascii');
xz=[pos -de];save([outfile '.xz'],'xz','-ascii');
%%
Poly.node=xz;
Poly.node(:,3)=0;
Poly.edge=(1:size(xz,1)-1)';Poly.edge(:,2)=Poly.edge(:,1)+1;
Poly.edge(:,3)=0;
elec=Out.elec;elec(:,3)=-99;
ln=size(Poly.node,1);
pbound=100;bound=2000;paradepth=15;
% pbound=50;bound=200; %test!
xmini=min(Poly.node(:,1))-pbound;
xmaxi=max(Poly.node(:,1))+pbound;
xmino=min(Poly.node(:,1))-bound;
xmaxo=max(Poly.node(:,1))+bound;
zmini=min(Poly.node(:,2))-paradepth;
zmino=min(Poly.node(:,2))-bound;
newnode=[xmaxi Poly.node(end,2);xmaxi zmini;xmini zmini;xmini Poly.node(1,2)];
newnode=[newnode;xmaxo 0;xmaxo zmino;xmino zmino;xmino 0];
newnode(:,3)=0;
Poly.node=[Poly.node;newnode;elec];
newedge=[0 1;1 2;2 3;3 4;4 1-ln;5 6;6 7;7 8;8 5]+ln;
newedge(:,3)=-1;newedge(6:8,3)=-2;
Poly.edge=[Poly.edge;newedge];
Poly.region=[xmino+1 zmino+1 1 0;xmini+1 zmini+1 2 0];
% figure(4);show2dpoly(Poly);axis normal;
writepoly2d([outfile '.poly'],Poly);
fid=fopen(fullfile(pname,'inv.cfg'),'w');
fprintf(fid,'DATAFILE=%s.ohm\n',fname);
fprintf(fid,'DIMENSION=2\n');
fprintf(fid,'SURFACESMOOTH=1\n');
fprintf(fid,'PARADX=0.2\n');
fprintf(fid,'PARA2DQUALITY=33.0\n');
fprintf(fid,'UNDERWATER=1\n');
fprintf(fid,'#PARAGEOMETRY=''cp %s.poly mesh/mesh.poly''\n',fname);
%fprintf(fid,'REFRAKTOR=%s.xz\n',fname);
fprintf(fid,'SPACECONFIG=2\n');
fprintf(fid,'RHOSTART=1.15\n');
fprintf(fid,'NOPROLONGATION=1\n');
fprintf(fid,'ZPOWER=0.5\n');
fprintf(fid,'INPUTERRVOLTAGE=4e-5\n');
fclose(fid);


% --- Executes on button press in inv1d.
function inv1d_Callback(hObject, eventdata, handles)
global Data DATA Geolore Log datafile Poly RES THK
set(gcf,'Pointer','watch');
null=datenum('00:00');
tmin=datenum(get(handles.editmintime,'String'),15)-null;
tmax=datenum(get(handles.editmaxtime,'String'),15)-null;
fi=find((Data.T>=tmin)&(Data.T<=tmax));
dx=sqrt(diff(DATA.lat(1:16:end)).^2+diff(DATA.lon(1:16:end)).^2);
DATA.x=[0;cumsum(dx)]*1000;
pos=interp1(DATA.T,DATA.x,Data.T(fi));
ch=1:16;ch(4)=[];
Sond=rmfield(delmeasurement(DATA.Sond,4),{'r','k','eind','names','nr'});
dl=50;% genauer heraus finden, u.U. abhängig von d selbst machen
vv=Geolore.v;vv(2:end+1)=(vv+Geolore.v)/2;vv(end)=Geolore.v(end);
dtime=dl./(vv/3.6);
% dep=interp1(Log.time,Log.depth,Geolore.gpstime-dtime);
% dep=interp1(Geolore.time,Geolore.depth,Geolore.gpstime-dtime);
dep=interp1(Geolore.time,Geolore.depth,Geolore.gpstime);
% edep=interp1(Geolore.x,dep,pos,'linear','extrap');
edep=interp1(Geolore.time,Geolore.depth,Data.T(fi),'linear','extrap');
RES=zeros(length(pos),6);THK=zeros(length(pos),5);
logdep=interp1(Log.time+mod(datenum('0:0:25'),1),Log.depth,Data.T(fi));
thk=[0 1 2 3 4];
for i=1:length(pos),
    Sond.elec(:,2)=edep(i);
    Sond.k=getkonf2d(Sond);
    Sond.rho=Data.R(i,ch)';
    Sond.r=Sond.rho.*Sond.k;
    writedcpfile(Sond,['calc\sond' num2str(i) '.dcp']);
    thk(1)=logdep(i);if isnan(logdep(i)), thk(1)=edep(i)+0.5; end
    writemodfile(i,'calc\sond.mod',[1.15],thk);
    cd('calc');system('em1dinv32 sond.mod >em1dinv.log');cd('..');
    res=emo2par('calc\sond.emo');
    RES(i,:)=res;THK(i,:)=thk;
end
figure(7);clf;
% cmin=log(min(RES(:)));cmax=log(max(RES(:)));
[cmin,cmax]=interperc(log(RES));
cmap=colormap(jet(64));
for i=1:size(RES,1),
    res=RES(i,:);thk=THK(i,:);
    d=[cumsum([0 thk]) 25];    
    for j=1:length(res), 
        cind=getcindex(log(res(j)),cmin,cmax,64);
        patch([-1 1 1 -1]*0.5+i,d(j+[0 0 1 1]),cmap(cind,:));
    end
end
hold on;plot(1:size(RES,1),edep,'w.');hold off;
set(gca,'YDir','reverse','XLim',[0.5 size(RES,1)+0.5]);
xt=get(gca,'XTick');if xt(1)>5, xt=[1 xt]; end
if xt(end)<size(RES,1)-5, xt=[xt size(RES,1)]; end
set(gca,'XTick',xt,'XTickLabel',datestr(Data.T(fi(xt)),15));
xlabel('time');ylabel('depth in m');
caxis([cmin,cmax]);cb=colorbar('horiz');
set(cb,'XTickLabel',num2strcell(rndig(exp(get(cb,'XTick')))));
set(gcf,'Pointer','arrow');