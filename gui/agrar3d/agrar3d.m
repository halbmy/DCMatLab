function varargout = agrar3d(varargin)
% AGRAR3D M-file for agrar3d.fig
%      AGRAR3D, by itself, creates a new AGRAR3D or raises the existing
%      singleton*.
%
%      H = AGRAR3D returns the handle to a new AGRAR3D or the handle to
%      the existing singleton*.
%
%      AGRAR3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AGRAR3D.M with the given input arguments.
%
%      AGRAR3D('Property','Value',...) creates a new AGRAR3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before agrar3d_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to agrar3d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help agrar3d

% Last Modified by GUIDE v2.5 15-Jan-2009 21:05:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @agrar3d_OpeningFcn, ...
                   'gui_OutputFcn',  @agrar3d_OutputFcn, ...
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


% --- Executes just before agrar3d is made visible.
function agrar3d_OpeningFcn(hObject, eventdata, handles, varargin)
global Sond
% Choose default command line output for agrar3d
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
A=load('default.conf');
Sond.ab2=A(:,1);
Sond.mn2=A(:,2);
Sond.c1=A(:,3);
Sond.ch=A(:,4);
Sond.k=pi./(1./abs(Sond.ab2-Sond.mn2)-1./abs(Sond.ab2+Sond.mn2));
Sond.sinn=find(A(:,end));

% UIWAIT makes agrar3d wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = agrar3d_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in readdatafile.
function readdatafile_Callback(hObject, eventdata, handles)
global datafile DATA
filetype='*.txt';infile=filetype;
if ~isempty(datafile),
    [pp,ff,ee]=fileparts(datafile);
    infile=fullfile(pp,infile);
end
[fname,pname]=uigetfile(filetype,'Load data file',infile);
if ~isstr(fname), return; end
datafile=fullfile(pname,fname); 
DATA=readlongresecsfile(datafile);
DATA.t=DATA.time(1:16:end);
fies=find(diff(DATA.t)<0);
if fies,
    st=sprintf('Sort data at %s, line %d',datestr(DATA.t(fies(1)),13),fies(1)*16);
    errordlg(st,'Data wrongly sorted');
    return;
    [DATA.t,II]=sort(DATA.t);
    DATA.lon=DATA.lon(II);
    DATA.lat=DATA.lat(II);
end
[DATA.y,DATA.x]=gps2xyz(DATA.lat(1:16:end),DATA.lon(1:16:end));
DATA.lon=(DATA.lon-floor(min(DATA.lon)))*1000;
DATA.lat=(DATA.lat-floor(min(DATA.lat)))*1000;
DATA.l=[0;cumsum(sqrt(diff(DATA.x).^2+diff(DATA.y).^2))];
set(handles.datastatus,'String',sprintf('Read %d single data from %s-%s',length(DATA.time),DATA.tstr{1},DATA.tstr{end}));
set(handles.tmin,'String',datestr(DATA.t(1),15));
set(handles.tmax,'String',datestr(DATA.t(end),15));
set(handles.figure1,'Name',['AGRAR3D - ' datafile]);


% --- Executes on button press in plotroute.
function plotroute_Callback(hObject, eventdata, handles)
global DATA
set(figure(10),'Name','Route','NumberTitle','Off');
% lo=DATA.lon(1:16:end);la=DATA.lat(1:16:end);
rr=24*60/str2num(get(handles.setmark,'String')); %round to 1 minute
mit=ceil(min(DATA.time)*rr);mat=floor(max(DATA.time)*rr);
tt=(mit:mat)/rr;
% t1=DATA.t;x1=DATA.x;
fi=find(diff(DATA.t)>0);
x1=interp1(DATA.t(fi),DATA.x(fi),tt);y1=interp1(DATA.t(fi),DATA.y(fi),tt);
x10=floor(mean(DATA.x)/1000)*1000;
y10=floor(mean(DATA.y)/1000)*1000;
x1=x1-x10;y1=y1-y10;
plot(DATA.x-x10,DATA.y-y10);axis xy equal tight;grid on;hold on;
for i=1:length(x1), 
    plot(x1(i),y1(i),'x');text(x1(i),y1(i),datestr(tt(i),15)); 
end;hold off
xlabel(sprintf('( x - %d )/m',round(x10)));
ylabel(sprintf('( y - %d )/m',round(y10)));
set(handles.routestatus,'String',sprintf('length = %.1fkm',DATA.l(end)/1000));

% --- Executes during object creation, after setting all properties.
function setmark_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function setmark_Callback(hObject, eventdata, handles)


% --- Executes on button press in loadconf.
function loadconf_Callback(hObject, eventdata, handles)
global Sond
[fname,pname]=uigetfile('*.conf','Load configuration');
confname=fullfile(pname,fname);
A=load(confname);
Sond.ab2=A(:,1);
Sond.mn2=A(:,2);
Sond.c1=A(:,3);
Sond.ch=A(:,4);
Sond.k=pi./(1./abs(Sond.ab2-Sond.mn2)-1./abs(Sond.ab2+Sond.mn2));
Sond.sinn=find(A(:,end));


% --- Executes on button press in rawdata.
function rawdata_Callback(hObject, eventdata, handles)
global DATA Sond
RHOA=ones(fix(length(DATA.u)/16)+1,length(Sond.sinn))*NaN;
ERR=RHOA;
for i=1:length(Sond.sinn),
    fi=find((DATA.ch==Sond.ch(Sond.sinn(i)))&(DATA.c1==Sond.c1(Sond.sinn(i))));
    R=DATA.u(fi)./DATA.i(fi);
    RHOA(fix(fi/16)+1,i)=DATA.u(fi)./DATA.i(fi)*Sond.k(Sond.sinn(i));
    ERR(fix(fi/16)+1,i)=DATA.d(fi);
end
RHOA(RHOA<=0)=NaN;RHOA(isinf(RHOA))=NaN;
% set(figure(1),'Name','Raw data','NumberTitle','Off','Renderer','painters');pause(0.1);
set(figure(1),'Name','Data','NumberTitle','Off','Renderer','painters');pause(0.1);
subplot(3,1,1);
LSIGMAA=log10(1000./RHOA);
imagesc(LSIGMAA');alpha(double(~isnan(RHOA))');
caxis(interperc(LSIGMAA));colormap(flipud(jet(64)));
set(gca,'XTickLabel',datestr(DATA.t(get(gca,'XTick')),15),'XTickMode','manual','XTickLabelMode','manual')
cb=colorbar;
set(cb,'YTickLabel',num2strcell(round(10.^get(cb,'YTick'))));
% set(cb,'XTick',mean(get(cb,'XLim')),'XTickLabel','\sigma_a in mS/m');
set(get(cb,'Title'),'String','\sigma_a in mS/m');
%%
% set(figure(2),'Name','Measured errors','NumberTitle','Off');pause(0.1);
subplot(3,1,2);
ERR(abs(ERR)>100)=NaN;
imagesc(log10(ERR'+1e-2));caxis(interperc(log10(ERR+1e-3)));
set(gca,'XTickLabel',datestr(DATA.t(get(gca,'XTick')),15),'XTickMode','manual','XTickLabelMode','manual')
cb=colorbar;
set(cb,'YTickLabel',num2strcell(round(10.^get(cb,'YTick'))));
set(get(cb,'Title'),'String','\epsilon in %');
% figure(1);
RHOA((ERR>10)|(RHOA>1100)|(RHOA<8))=NaN;
DATA.RHOA=RHOA(1:length(DATA.t),:);DATA.ERR=ERR(1:length(DATA.t),:);
% fi=find(diff(DATA.t)<=0)+1;
% DATA.x(fi)=[];DATA.y(fi)=[];DATA.t(fi)=[];


% --- Executes on button press in resample.
function resample_Callback(hObject, eventdata, handles)
global DATA Data
dx=str2num(get(handles.dx,'String'));
nc=str2num(get(handles.nc,'String'));
x10=round(mean(DATA.x)/1000)*1000;
y10=round(mean(DATA.y)/1000)*1000;
% v=diff(len)/diff(DATA.t);plot(DATA.t(1:end-1),v);datetick;
tmin=datenum(get(handles.tmin,'String'));
tmax=datenum(get(handles.tmax,'String'));
fmin=min(find(DATA.t>=tmin));
fmax=max(find(DATA.t<=tmax));

Data.l=(DATA.l(fmin):dx:DATA.l(fmax))';%newx
ff=[1;find(diff(DATA.l)>0)+1];
Data.t=interp1(DATA.l(ff),DATA.t(ff),Data.l); %newt
fi=find(diff(DATA.t)>0);
Data.x=interp1(DATA.t(fi),DATA.x(fi),Data.t,'linear','extrap');
Data.y=interp1(DATA.t(fi),DATA.y(fi),Data.t,'linear','extrap');
Data.x0=floor(mean(Data.x)/1000)*1000;
Data.y0=floor(mean(Data.y)/1000)*1000;
Data.x=Data.x-Data.x0;
Data.y=Data.y-Data.y0;
%%
rhol=1;rhou=2000; %in GUI integrieren!!!
Data.RHOA=zeros(length(Data.t),size(DATA.RHOA,2));
RHOA=DATA.RHOA(fmin:fmax,:);ERR=DATA.ERR(fmin:fmax,:);DL=DATA.l(fmin:fmax);
for i=1:size(DATA.RHOA,2),
    fi=find(isfinite(RHOA(:,i)));
    datain=logtrans(RHOA(fi,i),rhol,rhou);
    dataout=harmfit(DL(fi),datain,nc,Data.l,ERR(fi,i)/100+0.3);
    if i==0,%6
        figure(22);
        d1=(ERR(fi,i)/100+1).*datain;
        d2=(-ERR(fi,i)/100+1).*datain;
        ll=DL(fi);
        plot(ll,datain,'b-',Data.l,dataout,'r-');%,ll,d1,'g-',ll,d2,'g-',);
    end
    Data.RHOA(:,i)=ilogtrans(dataout,rhol,rhou);
end
% set(figure(3),'Name','Processed data','NumberTitle','Off');
figure(1);subplot(3,1,3);
imagesc(Data.l,1:size(DATA.RHOA,2),log10(1000./Data.RHOA)');
colormap(flipud(jet(64)));
% datetick;
axis tight;cb=colorbar;
set(cb,'YTickLabel',num2strcell(round(10.^get(cb,'YTick'))));
set(get(cb,'Title'),'String','\sigma_a in mS/m');


% --- Executes during object creation, after setting all properties.
function dx_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function dx_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function nc_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function nc_Callback(hObject, eventdata, handles)


% --- Executes on button press in inversion.
function inversion_Callback(hObject, eventdata, handles)
global Data Model Sond
oldpwd=pwd;
try,
Model.x=Data.x;Model.y=Data.y;
Model.x0=Data.x0;Model.y0=Data.y0;
Model.RES=ones(size(Data.RHOA,1),5);
DCP=[Sond.mn2(Sond.sinn) Sond.ab2(Sond.sinn) Data.RHOA(1,:)'];
nsond=size(Data.RHOA,1);
wb=-1;
if 0,
    DD=abmn2n(-Sond.ab2(Sond.sinn),Sond.ab2(Sond.sinn),-Sond.mn2(Sond.sinn),Sond.mn2(Sond.sinn));
    saveinv2dfile('sond.dat',DD);
    ss='';for i=1:size(Data.RHOA,2), ss=[ss '%f\t']; end
    fid=fopen('sand.rhoa','w');
    ss(end)='n';fprintf(fid,ss,Data.RHOA');
    fclose(fid);
    display('dc1dlong -v -e5 -d sond.dat sand.rhoa');    
    Model.RES=loadsens('model.matrix');
else
    nupd=floor(nsond/25)+1;
    lupd=0;
    wb=waitbar(0,'1D Inversion');
    tic;
    for i=1:nsond,
        DCP(:,3)=Data.RHOA(i,:)';
        % writedcpfile(DCP,['calc\sond' num2str(i) '.dcp']);
        writedcpfile(DCP,'calc\sond.dcp');
        cd('calc');system('em1dinv32 sond.mod >em1dinv.log');cd('..');
        res=emo2par('calc\sond.emo');
        Model.RES(i,:)=res;
        lupd=lupd+1;
        if lupd>=nupd,
            waitbar(i/nsond,wb); 
            lupd=0;
        end
    end
    close(wb);
    toc;
end
catch,
    cd(oldpwd);
    if ishandle(wb), close(wb); end    
end
agrar3d('plotresults_Callback',gcbo,[],guidata(gcbo)); 

% --- Executes on button press in plotresults.
function plotresults_Callback(hObject, eventdata, handles)
global Model Data
set(figure(4),'Name','Inversion Results','NumberTitle','Off');
%% plot results
dx=5;
dx=str2num(get(handles.dx,'String'));
xi=min(Model.x):dx:max(Model.x);yi=min(Model.y):dx:max(Model.y);
[XI,YI]=meshgrid(xi,yi);cax=fliplr(3-interperc(log10(Model.RES)));
warning off MATLAB:griddata:DuplicateDataPoints
for i=1:5,
    RI=10.^griddata(Model.x,Model.y,log10(Model.RES(:,i)),XI,YI);
%     set(figure(10+i),'Name',['Layer' num2str(i)],'NumberTitle','Off');
    subplot(2,3,i);
    imagesc(xi,yi,log10(1000./RI));alpha(double(~isnan(RI)));
    colormap(flipud(jet(16)));
    caxis(cax);%cb=colorbar;
%     set(cb,'YTickLabel',num2strcell(round(1000./10.^get(cb,'YTick'))),'XTick',0.5,'XTickLabel','EC[mS/m]');
    hold on;plot(Model.x,Model.y,'kx-','MarkerSize',1);hold off;axis xy equal tight;
end
subplot(2,3,6);
cbar(10.^cax(1),10.^cax(2),1,1);
title('\sigma_a in mS/m');


% --- Executes on button press in exportfigures.
function exportfigures_Callback(hObject, eventdata, handles)
global datafile modelfile
if ishandle(10), epsprint(10,strrep(datafile,'.txt','-route'),1); end
% if ishandle(1), exportpng(1,strrep(datafile,'.txt','-rawdata.png')); end
if ishandle(1), exportpng(1,strrep(datafile,'.txt','-data.png')); end
if ishandle(2), epsprint(2,strrep(datafile,'.txt','-errors'),1); end
if ishandle(3), epsprint(3,strrep(datafile,'.txt','-filtdata'),1); end
if ishandle(4), 
    figure(4);pause(0.3);
    if ~isempty(modelfile), exfile=strrep(modelfile,'.model',''); 
    else exfile=strrep(datafile,'.txt',''); end
    epsprint(4,[exfile '-result'],1);
    exportpng(4,[exfile '-result.png']);
end


% --- Executes on button press in savemodel.
function savemodel_Callback(hObject, eventdata, handles)
global Model Data datafile modelfile
model=[Model.x+Model.x0 Model.y+Model.y0 1000./Model.RES]';
outfile=strrep(datafile,'.txt','.model');
[fname,pname]=uiputfile('*.model','Model file',outfile);
if isstr(fname),
    [pp,nn,ff]=fileparts(fname);
    if isempty(ff), fname=[fname '.model']; end
    modelfile=fullfile(pname,fname);
    fid=fopen(modelfile,'w');
    fprintf(fid,'x[m]\ty[m]\tEC1[mS/m]\tEC2[mS/m]\tEC3[mS/m]\tEC4[mS/m]\tEC5[mS/m]\r\n');
    fprintf(fid,'%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\r\n',model);
    fclose(fid);
end


% --- Executes on button press in loadmodel.
function loadmodel_Callback(hObject, eventdata, handles)
global Model datafile modelfile
filetype='*.model';infile=filetype;
if ~isempty(datafile),
    [pp,ff,ee]=fileparts(datafile);
    infile=fullfile(pp,infile);
end
[fname,pname]=uigetfile(filetype,'Load data file',infile);
if ~isstr(fname), return; end
modelfile=fullfile(pname,fname); 
A=textread(modelfile,'','headerlines',1);
Model.x=A(:,1);Model.y=A(:,2);
Model.RES=1000./A(:,3:end);
Model.x0=floor(mean(Model.x)/1000)*1000;
Model.y0=floor(mean(Model.y)/1000)*1000;
Model.x=Model.x-Model.x0;
Model.y=Model.y-Model.y0;
datafile=strrep(modelfile,'.model','.txt');
agrar3d('plotresults_Callback',gcbo,[],guidata(gcbo)); 

% --- Executes on button press in addmodel.
function addmodel_Callback(hObject, eventdata, handles)
global Model datafile
filetype='*.model';infile=filetype;
if ~isempty(datafile),
    [pp,ff,ee]=fileparts(datafile);
    infile=fullfile(pp,infile);
end
[fname,pname]=uigetfile(filetype,'Load data file',infile);
if ~isstr(fname), return; end
modelfile=fullfile(pname,fname); 
A=textread(modelfile,'','headerlines',1);
Model.x=[Model.x;A(:,1)-Model.x0];
Model.y=[Model.y;A(:,2)-Model.y0];
Model.RES=[Model.RES;1000./A(:,3:end)];
agrar3d('plotresults_Callback',gcbo,[],guidata(gcbo)); 


% --- Executes during object creation, after setting all properties.
function tmin_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function tmin_Callback(hObject, eventdata, handles)
global DATA
tmin=datenum(get(handles.tmin,'String'));
tmax=datenum(get(handles.tmax,'String'));
fmin=min(find(DATA.t>=tmin));
fmax=max(find(DATA.t<=tmax));
if ishandle(1), 
    figure(1);
%     set(gca,'Xlim',[fmin fmax]); 
    ch=get(gcf,'Children');
    for i=1:length(ch), 
        xl=get(ch(i),'Xlim');
        if xl(1)>0, set(ch(i),'Xlim',[fmin fmax]); end
    end
end

% --- Executes during object creation, after setting all properties.
function tmax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function tmax_Callback(hObject, eventdata, handles)
global DATA
tmin=datenum(get(handles.tmin,'String'));
tmax=datenum(get(handles.tmax,'String'));
fmin=min(find(DATA.t>=tmin));
fmax=max(find(DATA.t<=tmax));
if ishandle(1), 
    figure(1);%subplot(3,1,1);
%     set(gca,'Xlim',[fmin fmax]); 
    ch=get(gcf,'Children');
    for i=1:length(ch), 
        xl=get(ch(i),'Xlim');
        if xl(1)>0, set(ch(i),'Xlim',[fmin fmax]); end
    end
end


