function varargout = underwater2(varargin)
% UNDERWATER2 M-file for underwater2.fig
%      UNDERWATER2, by itself, creates a new UNDERWATER2 or raises the existing
%      singleton*.
%
%      H = UNDERWATER2 returns the handle to a new UNDERWATER2 or the handle to
%      the existing singleton*.
%
%      UNDERWATER2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UNDERWATER2.M with the given input arguments.
%
%      UNDERWATER2('Property','Value',...) creates a new UNDERWATER2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before underwater2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to underwater2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help underwater2

% Last Modified by GUIDE v2.5 22-Feb-2008 14:39:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @underwater2_OpeningFcn, ...
                   'gui_OutputFcn',  @underwater2_OutputFcn, ...
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


% --- Executes just before underwater2 is made visible.
function underwater2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to underwater2 (see VARARGIN)

% Choose default command line output for underwater2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes underwater2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = underwater2_OutputFcn(hObject, eventdata, handles) 
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
if i==7,
    [All{1},All{2},All{3},All{4},All{5},All{6},All{7}]=textread(datafile,formstr,'headerlines',1);
elseif i==8,
    [All{1},All{2},All{3},All{4},All{5},All{6},All{7},All{8}]=textread(datafile,formstr,'headerlines',1);
else
    display('Number of fields does not match!');return;
end
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
    t0=datenum('02:00:00');
    DATA.time=datenum(DATA.tstr)+t0;
    DATA.time=DATA.time-floor(median(DATA.time));
end
set(handles.statusdatafile,'String',sprintf('%d data (%s-%s)',...
    DATA.ndata,datestr(min(DATA.time),15),datestr(max(DATA.time),15)));
adcfile=strrep(datafile,'.txt','.ADC');
if exist(adcfile,'file'),
    set(handles.editgeolorefile,'String',adcfile);
    underwater2('editgeolorefile_Callback',gcbo,[],guidata(gcbo));
end    
logfile=strrep(datafile,'.txt','.log');
if exist(logfile,'file'),
    set(handles.editdepthlogfile,'String',logfile);
    underwater2('editdepthlogfile_Callback',gcbo,[],guidata(gcbo));
end    
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
underwater2('editdatafile_Callback',gcbo,[],guidata(gcbo));


function editgeolorefile_Callback(hObject, eventdata, handles)
global Geolore
geolorefile=get(handles.editgeolorefile,'String');
Geolore=[];
try,
fid=fopen(geolorefile);
zeile=fgetl(fid);
fclose(fid);
if strfind(zeile,',N'), %GPS data valid
    [Time,D,T,Geolore.LF,Geolore.lat,Geolore.lon]=textread(geolorefile,'%*s%s%f%f%f%f,N%f,E%*s%*s');
else %missing
%     [Time,D,T,Geolore.LF,Geolore.lat,Geolore.lon]=textread(geolorefile,'%*s%s%f%f%f%f,*s%f,*s%*s%*s');
    [Time,D,T,Geolore.LF]=textread(geolorefile,'%*s%s%f%f%f , ,%*s%*s');
end
Geolore.time=datenum(Time)+datenum('2:0:0');
Geolore.time=Geolore.time-floor(median(Geolore.time));
Geolore.depth=D-min(min(D),0);
if median(diff(Geolore.time))<4/86400,
    Geolore.time=Geolore.time(1:5:end);
    Geolore.depth=Geolore.depth(1:5:end);
    if isfield(Geolore,'lat'), Geolore.lat=Geolore.lat(1:5:end); end
    if isfield(Geolore,'lon'), Geolore.lon=Geolore.lon(1:5:end); end
    Geolore.LF=Geolore.LF(1:5:end);
end
if isfield(Geolore,'lat')&isfield(Geolore,'lon'),
    [Geolore.y,Geolore.x]=gps2xyz(Geolore.lat,Geolore.lon);
    % fl=10.^floor(log10(max([diff(minmax(Geolore.x)) diff(minmax(Geolore.x))])));
    % Geolore.x=(Geolore.x/fl-floor(mean(Geolore.x)/fl))*fl;
    % Geolore.y=(Geolore.y/fl-floor(mean(Geolore.y)/fl))*fl;
    dl=sqrt(diff(Geolore.x).^2+diff(Geolore.y).^2);
    Geolore.v=dl./diff(Geolore.time)/24/1000/1.852; %kn
    Geolore.l=cumsum([0;dl]);
else
    Geolore.v=4.5;
end
set(handles.statusgeolorefile,'String',sprintf('%d readings (%s-%s), d=%.1f-%.1f',...
    length(Geolore.time),datestr(min(Geolore.time),15),datestr(max(Geolore.time),15),min(Geolore.depth),max(Geolore.depth)));
catch,
    uiwait(errordlg({'Error reading data logger file',lasterr},'Error'));
end

% --- Executes during object creation, after setting all properties.
function editgeolorefile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in browsegeolorefile.
function browsegeolorefile_Callback(hObject, eventdata, handles)
global datafile
filetype='*.ADC';infile=filetype;
if ~isempty(datafile),
    [pp,ff,ee]=fileparts(datafile);
    infile=fullfile(pp,infile);
end
[fname,pname]=uigetfile(filetype,'Load Geolore file',infile);
if ~isstr(fname), return; end
geolorefile=fullfile(pname,fname);
set(handles.editgeolorefile,'String',geolorefile);
underwater2('editgeolorefile_Callback',gcbo,[],guidata(gcbo));

function editdepthlogfile_Callback(hObject, eventdata, handles)
global Log
depthlogfile=get(handles.editdepthlogfile,'String');
% [HH,MM,SS,HU,D15,D100]=textread(depthlogfile,'%d:%d:%d.%d,%d,%d\r\n%*c');
[HH,MM,SS,HU,D15,D100]=textread(depthlogfile,'%d:%d:%d.%d,%d,%d\r\n');
Log.time=(((SS+HU/100)/60+MM)/60+HH)/24;
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
filetype='*.log';infile=filetype;
if ~isempty(datafile),
    [pp,ff,ee]=fileparts(datafile);
    infile=fullfile(pp,infile);
end
[fname,pname]=uigetfile(filetype,'Load depth log file',infile);
if ~isstr(fname), return; end
gpslogfile=fullfile(pname,fname);
set(handles.editdepthlogfile,'String',gpslogfile);
underwater2('editdepthlogfile_Callback',gcbo,[],guidata(gcbo));



% --- Executes on button press in showroute.
function showroute_Callback(hObject, eventdata, handles)
global DATA Geolore
figure(11);
%%
% plot(DATA.lon,DATA.lat,'b.');
[y,x]=gps2xyz(DATA.lat,DATA.lon);
plot(x,y,'b.');
axis equal
grid on
if isfield(Geolore,'lon'),
%     hold on;plot(Geolore.lon,Geolore.lat,'r.');hold off;
    hold on;plot(Geolore.x,Geolore.y,'r.');hold off;
end
tt=linspace(min(DATA.time),max(DATA.time),20);
% tlon=interp1(DATA.time(1:16:end),DATA.lon(1:16:end),tt);
% tlat=interp1(DATA.time(1:16:end),DATA.lat(1:16:end),tt);
% for i=1:length(tt),    
%     set(text(tlon(i),tlat(i),datestr(tt(i),15)),'HorizontalAlignment','left');
% end
tx=interp1(DATA.time(1:16:end),x(1:16:end),tt);
ty=interp1(DATA.time(1:16:end),y(1:16:end),tt);
for i=1:length(tt),    
    set(text(tx(i),ty(i),datestr(tt(i),15)),'HorizontalAlignment','left');
end
legend('data','logger');
%%

% --- Executes on button press in showdepths.
function showdepths_Callback(hObject, eventdata, handles)
global Geolore Log
figure(12);clf;
%%
mul=1;%mul=1/1.852;
if isfield(Geolore,'v'),
   subplot(2,1,1);
   plot(Geolore.time(1:end-1)+diff(Geolore.time)/2,Geolore.v*mul,'b-');%
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
axis ij;xlabel('time');ylabel('z in m');datetick;grid on;
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
DATA.Sond=readsnd2d('sond2.s2d');
DATA.Sond.elec(:,2)=8;
if isfield(Geolore,'depth'), DATA.Sond.elec(:,2)=median(Geolore.depth); end
DATA.Sond.k=getkonf2d(DATA.Sond);
DATA.Rhoa0=DATA.R*diag(DATA.Sond.k);
DATA.T=DATA.time(1:16:end);
if length(DATA.T)>size(DATA.R,1), DATA.T=DATA.T(1:size(DATA.R,1)); end
DATA.Rhoa0(DATA.Rhoa0<0)=NaN;
DATA.Rhoa0(DATA.Rhoa0>interperc(DATA.Rhoa0,99.9))=NaN;
DATA.Rhoa0(DATA.Err>20)=NaN;
DATA.Rhoa0(:,6)=NaN;
DATA.R(isnan(DATA.Rhoa0))=NaN;
DATA.cax=10.^interperc(log10(DATA.Rhoa0),[2 98]);
DATA.cax=[1 3];
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
%% dt einlesen
dt=10/86400;
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
Data.Rhoa(:,6)=NaN;
Data.Rhoa(Data.Rhoa<=0)=NaN;
figure(3);imageline(Data.T,Data.Rhoa,Data.cax);
set(handles.editmintime,'String',datestr(min(Data.T),15));
set(handles.editmaxtime,'String',datestr(max(Data.T),15));


% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)
global datafile
names={'raw','rhoa0','rhoa','','','','1dinv','','','','route','depth'};
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
pos=interp1(Geolore.time,Geolore.l,Data.T(fi));
%%
ch=1:16;ch(6)=[];
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
% dep=interp1(Geolore.time,Geolore.depth,Geolore.gpstime);
dep=Geolore.depth;
[unl,fifi]=unique(Geolore.l);
Out.elec(:,2)=interp1(unl,dep(fifi),Out.elec(:,1),'linear','extrap');
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
dtime=datenum('0:0:25');dtime=dtime-floor(dtime);
dtime=datenum('0:0:40');dtime=dtime-floor(dtime);
de=interp1(Log.time+dtime,Log.depth,Data.T(fi),'linear','extrap');
% plot(pos,de,'r-',Out.elec(:,1),Out.elec(:,2),'b-');axis ij;legend('log','sensor');
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
fprintf(fid,'UNDERWATER2=1\n');
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
if length(DATA.x)>length(DATA.T), DATA.x=DATA.x(1:length(DATA.T)); end
pos=interp1(DATA.T,DATA.x,Data.T(fi));
ch=1:16;ch(6)=[];
Sond=rmfield(delmeasurement(DATA.Sond,6),{'r','k','eind','names','nr'});
dl=25;% genauer heraus finden, u.U. abhängig von d selbst machen
vv=Geolore.v;vv(2:end+1)=(vv+Geolore.v)/2;vv(end)=Geolore.v(end);
dtime=dl./(vv/3.6*1.852)/86400/2; %Zeitdifferenz zw. Datenlogger und Mittelpunkt
% dep=interp1(Log.time,Log.depth,Geolore.gpstime-dtime);
% dep=interp1(Geolore.time,Geolore.depth,Geolore.gpstime-dtime);
% dep=interp1(Geolore.time,Geolore.depth,Geolore.gpstime);
dep=Geolore.depth;
% edep=interp1(Geolore.l,dep,pos,'linear','extrap');
edep=interp1(Geolore.time+dtime,Geolore.depth,Data.T(fi),'linear','extrap');
RES=zeros(length(pos),6);THK=zeros(length(pos),5);
% logdep=interp1(Log.time+mod(datenum('0:0:25'),1),Log.depth,Data.T(fi));
dl2=55;
% eigentlich vv auf Log.time interpolieren!!!
etime=dl2/(median(vv)/3.6*1.852)/86400/2;
datestr(median(dtime),13)
datestr(etime,13)
logdep=interp1(Log.time+etime,Log.depth,Data.T(fi));
% thk=[0 1 2 3 4];
thk=[0 0.5 1 2 3];
lf=median(Geolore.LF);
lfs=interp1(Geolore.time,Geolore.LF,Data.T(fi));
for i=1:length(pos),
    Sond.elec(:,2)=edep(i);
    Sond.k=getkonf2d(Sond);
    Sond.rho=Data.R(i,ch)';
    Sond.r=Sond.rho.*Sond.k;
    writedcpfile(Sond,['calc\sond' num2str(i) '.dcp']);
    thk(1)=logdep(i);if isnan(logdep(i)), thk(1)=edep(i)+0.5; end
    writemodfile(i,'calc\sond.mod',10/lfs(i),thk);
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
    d=[cumsum([0 thk])];d(end+1)=d(end)+thk(end);
    for j=1:length(res), 
        cind=getcindex(log(res(j)),cmin,cmax,64);
        patch([-1 1 1 -1]*0.5+i,d(j+[0 0 1 1]),cmap(cind,:),'LineStyle','none');
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
set(cb,'YTick',mean(get(cb,'Ylim')),'YTickLabel','Ohmm');
set(gcf,'Pointer','arrow');