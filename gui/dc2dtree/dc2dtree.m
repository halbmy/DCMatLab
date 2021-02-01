function varargout = dc2dtree(varargin)

% Last Modified by GUIDE v2.5 12-Apr-2007 14:14:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dc2dtree_OpeningFcn, ...
                   'gui_OutputFcn',  @dc2dtree_OutputFcn, ...
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


% --- Executes just before dc2dtree is made visible.
function dc2dtree_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for dc2dtree
handles.output = hObject;

global libmmfile
libmmfile=checklict;
if ~isequal(libmmfile,4), 
    uiwait(errordlg('No valid license file found!'));
end
% Update handles structure
guidata(hObject, handles);
iconify(hObject);

% UIWAIT makes dc2dtree wait for user response (see UIRESUME)
% uiwait(handles.figure1);
global ana
%ana=load('48.pot','-ascii');
fid=fopen('48.pot','r');
ana=fscanf(fid,'%f');
fclose(fid);
global progdir
progdir=pwd;
lam=30;
set(handles.lamtext,'String',num2str(lam));
set(handles.lamslide,'Val',1-log10(lam)/3);
set(handles.figure1,'Name','DC2dTree  v.0.9.2');
set(handles.message,'String',{'DC2dTree - Tree Tomography';'version 0.9.2';'------------RESISTIVITY.NET------------'});
set(handles.run,'Enable','off');
set(handles.loaddata,'Enable','off');
set(handles.preview,'Enable','off');
set(handles.prepro,'Enable','off');
set(handles.optinv,'Enable','on');%off
if 0&&exist('start.mat','file'),
    load('start.mat');
    axes(handles.plot);
    tripatchmod(Mesh,res,N);
end
set(handles.plot,'Color',get(handles.figure1,'Color'),'XTick',[],'YTick',[]);
if 0&&exist('mehl150.mat','file'),
    load('mehl150.mat');
    axes(handles.plot);
    image(A);
    axis equal tight
    set(handles.plot,'XTick',[],'YTick',[]);
end

% --- Outputs from this function are returned to the command line.
function varargout = dc2dtree_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
% end

% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
dc2dtree('loadpro_Callback',gcbo,[],guidata(gcbo))
% end

% --- Executes on button press in prepro.
function prepro_Callback(hObject, eventdata, handles)
global Pro Mesh N progdir S Meshopts
set(handles.run,'Enable','Off');
set(handles.preview,'Enable','Off');
Pro.lambda=str2double(get(handles.lamtext,'String'));
if isfield(Pro,'dirname'), cd(Pro.dirname); end
try,
cmdline=[progdir filesep 'dc2dtreepre -v '];
if ~isfield(Meshopts,'nospline')||(Meshopts.nospline==0), cmdline=[cmdline '-B ']; end
if isfield(Meshopts,'equirefine')&&(Meshopts.equirefine>0), cmdline=[cmdline '-E ']; end
if isfield(Meshopts,'parasmooth')&&(Meshopts.parasmooth>0), cmdline=[cmdline '-S ']; end
if isfield(Meshopts,'paraquality'), cmdline=[cmdline '-Q' num2str(Meshopts.paraquality) ' ']; end
if isfield(Meshopts,'pararefine'), cmdline=[cmdline '-R' num2str(Meshopts.pararefine) ' ']; end
if isfield(Meshopts,'paramaxarea'), cmdline=[cmdline '-A' num2str(Meshopts.paramaxarea) ' ']; end
if isfield(Meshopts,'primsmooth')&&(Meshopts.primsmooth>0), cmdline=[cmdline '-s ']; end
if isfield(Meshopts,'primquality'), cmdline=[cmdline '-q' num2str(Meshopts.primquality) ' ']; end
if isfield(Meshopts,'primrefine'), cmdline=[cmdline '-r' num2str(Meshopts.primrefine) ' ']; end
if isfield(Meshopts,'primmaxarea'), cmdline=[cmdline '-a' num2str(Meshopts.primmaxarea) ' ']; end
if isfield(Meshopts,'secmeshrefine'), cmdline=[cmdline '-f' num2str(Meshopts.secmeshrefine) ' ']; end
cmdline=[cmdline '"' Pro.datfile '"'];
set(gcf,'Pointer','watch');
systemcall(cmdline);
set(gcf,'Pointer','arrow');
[fpath,fname,fext]=fileparts(Pro.datfile);
datafile=[fname '.data'];
% datafile=strrep(Pro.datfile,'.ohm','') '.data'];
if exist(datafile,'file'), 
    Pro.datafile=datafile; 
    appendmessage(handles,['Found corrected data file ' Pro.datafile]);
%     NN=readinv2dfile(Pro.datafile,1);N.rhoa=NN.r;N.k=NN.err*100;
%     N=readinv2dfile(Pro.datafile,1);
else
    appendmessage(handles,'Data file not found!');cd(progdir); return;
end
if exist('tmp\meshPara.bms','file'),%||exist('tmp\meshPara.e','file')&exist('tmp\meshPara.n','file'),
    Mesh=loadmesh('tmp\meshPara.bms');
    ss=sprintf('Found mesh %d elements, %d nodes',Mesh.ncells,Mesh.nnodes);
    appendmessage(handles,ss);
    Pro.res=ones(Mesh.ncells,1)*median(N.r);%rhoa
else
    appendmessage(handles,'No mesh found!');return;
end
if ~exist('tmp\sensMat.mat','file'),
    appendmessage(handles,'No sensitivity matrix found!');return; end
S=loadsens;
if ~exist('tmp\secPot.mat','file'),
    appendmessage(handles,'No potential matrix found!');return; end
if ~exist('tmp\secPot.collect','file'),
    appendmessage(handles,'No collect found!');return; end
cd(progdir);
if isfield(Pro,'datafile'),
    fdf=fullfile(Pro.dirname,Pro.datafile);
    if exist(fdf,'file'), N=readinv2dfile(fdf,1); end
end
if (~isfield(N,'err'))||(min(N.err)<=0), N.err=ones(size(N.a))*0.03; end
set(handles.run,'Enable','On');
set(handles.preview,'Enable','On');
catch,
    disp('Catched in prepare:');
    disp(lasterr);
    cd(progdir);
end
dc2dtree('savepro_Callback',gcbo,[],guidata(gcbo));
dc2dtree('showmod_Callback',gcbo,[],guidata(gcbo))
for i=1:size(N.elec,1), set(text(N.elec(i,1),N.elec(i,2),num2str(i)),...
        'Color','red','VerticalAlignment','middle','HorizontalAlignment','center'); end
% end

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
global Pro progdir Invopts
if isfield(Pro,'res'),
    Pro.res(:)=1;
    dc2dtree('showmod_Callback',gcbo,[],guidata(gcbo))
    drawnow;
end
set(handles.preview,'Value',0);
Pro.lambda=str2num(get(handles.lamtext,'String'));
cd(Pro.dirname);
try,
if ~isfield(Pro,'datafile')||~exist(Pro.datafile,'file'),
    uiwait(errordlg('No valid data file specified!'));cd(progdir);return;
end
delete('tmp\model_iter.*');
delete('tmp\modelResponse.*');
delete('tmp\modelReciprocity.*');
cmdline=[progdir filesep 'dc2dtreerun -v '];
if get(handles.optlambda,'Value')>0,
    cmdline=[cmdline '-O '];
else
    cmdline=[cmdline '-l' num2str(Pro.lambda) ' '];
end
if isfield(Invopts,'lowerbound'), cmdline=[cmdline '-b' num2str(Invopts.lowerbound) ' ']; end
if isfield(Invopts,'upperbound'), cmdline=[cmdline '-y' num2str(Invopts.upperbound) ' ']; end
if get(handles.robust,'Value')>0, cmdline=[cmdline '-R ']; end
if get(handles.blockymodel,'Value')>0, cmdline=[cmdline '-B ']; end
cmdline=[cmdline '"' Pro.datafile '"'];
set(gcf,'Pointer','watch');
systemcall(cmdline);
set(gcf,'Pointer','arrow');
if exist('tmp\chi2.vector','file'),
    fid=fopen('tmp\chi2.vector');chiq=fscanf(fid,'%f');fclose(fid);
    chiq=chiq(find(chiq));
    % appendmessage(handles,['CHI^2 fit: ' sprintf('-->%.1f',chiq)]);
    appendmessage(handles,sprintf('CHI^2 fit %.1f-->(%d It.)-->%.1f',chiq(1),length(chiq)-1,chiq(end)));
end
Pro.resfile=[strrep(Pro.name,'+','_') '.res'];
aa=dir('tmp\model_iter.*');
if isempty(aa),
    appendmessage(handles,'No model vector found!');
    return; 
end
for i=1:length(aa), dd{i}=aa(i).date; end
cd(progdir);
idx=sortcellchar(dd);
cd(Pro.dirname);
fprintf(['taking model file ' aa(idx(end)).name]);
if ispc,
  dos(['copy /y "tmp' filesep aa(idx(end)).name '" "' Pro.resfile '"']);
else
  dos(['cp "tmp' filesep aa(idx(end)).name '" "' Pro.resfile '"']);
end
Pro.morfile=[Pro.name '.mor'];
aa=dir('tmp\modelResponse.*');
if ispc,
  dos(['copy /y "tmp' filesep aa(idx(end)).name '" "' Pro.morfile '"']);
else
  dos(['copy /y "tmp' filesep aa(idx(end)).name '" "' Pro.morfile '"']);
end
fid=fopen(Pro.resfile);Pro.res=fscanf(fid,'%f');fclose(fid);
fid=fopen(Pro.morfile);Pro.mor=fscanf(fid,'%f');fclose(fid);
cd(progdir);
dc2dtree('showmod_Callback',gcbo,[],guidata(gcbo))
dc2dtree('savepro_Callback',gcbo,[],guidata(gcbo));
catch,
    disp('Catched in run procedure');
    disp(lasterr);
    cd(progdir);
end
% end

% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function newpro_Callback(hObject, eventdata, handles)
global N Pro progdir malstat ana S Mesh
set(handles.message,'String',{'DC2dTree - Tree Tomography';'-----------------------------------'});
alltypes='*.ddz;*.ohm;*.dat;*.tx0';
if isfield(Pro,'dirname')&&~isempty(Pro.dirname), 
    ss=strfind(Pro.dirname,filesep);
    infile=fullfile(Pro.dirname(1:ss(end)-1),alltypes);
else
    infile=alltypes;
end
if isfield(Pro,'dirname'), infile=[Pro.dirname filesep '..' filesep alltypes]; end
[fname,pname]=uigetfile({alltypes,'Known files';
    '*.ddz','Geotom Files(*.ddz)';'*.ohm','OHM Files(*.ohm)';...
    '*.dat','DAT Files(*.dat)';'*.tx0','Lippmann Files(*.tx0)';...
    '*.*','All files(*.*)'},'Load data file',infile);
if isequal(fname,0), return; end
set(handles.prepro,'Enable','Off');
set(handles.run,'Enable','Off');
set(handles.preview,'Enable','Off');
set(handles.robust,'Value',0.0);
set(handles.blockymodel,'Value',0.0);
dc2dtree('lamdefault_Callback',gcbo,[],guidata(gcbo))
Pro=[];
[ffp,ffn,ffe]=fileparts(fname);
switch ffe,
    case '.ddz',
        Pro.ddzfile=fullfile(pname,fname);
        appendmessage(handles,['Opening ddz file ' Pro.ddzfile]);
        ss=questdlg('Specify direction of electrode arrangement','Direction','Clockwise','Counterclockwise','Counterclockwise');
        if isequal(ss,'Counterclockwise'), di=-1; else di=1; end
        [N,Pro.radius]=loadddz(Pro.ddzfile,di);
        N.k=kfactorana(N,ana);
        appendmessage(handles,sprintf('Tree radius %.2m',Pro.radius));
    case '.tx0',
        Pro.ddzfile=fullfile(pname,fname);
        appendmessage(handles,['Opening tx0 file ' Pro.ddzfile]);
        [N,Pro.freq]=readtx0file(Pro.ddzfile);
        [fpit,ppit]=uigetfile('*.pit','Specify geometry file',fullfile(pname,'*.pit'));
        Shot=readpit(fullfile(ppit,fpit));
        if min(N.err)<0.01, N.err=N.err+0.02; end
        N.elec=Shot.pos;
        N.k=kfactorana(N,ana);
        N.r=abs(N.rho.*N.k);
    case '.dat',
        Pro.datfile=fullfile(pname,fname);
        N=sort2delecs(readinv2dfile(Pro.datfile,1));
        N.k=kfactorana(N,ana);
        if ~isfield(N,'rho'), N.rho=N.r; end
        N.r=abs(N.rho.*N.k);
    case '.ohm',
        N=sort2delecs(readinv2dfile(fullfile(pname,fname),1));
        Pro.datfile=fname;
        N.k=kfactorana(N,ana);
        if ~isfield(N,'rho'), N.rho=N.r; end
        N.r=abs(N.rho.*N.k);
    otherwise,
        display('Could not read data');
end
appendmessage(handles,sprintf('Loaded %d data with %d electrodes',length(N.a),size(N.elec,1)));
axes(handles.plot);
showpseudo(N,N.r);malstat=2;
if isequal(ffe,'.ddz'),
    answer=questdlg(['Radius is ' num2str(Pro.radius)],'Radius',...
        'Keep','Change','Load File','Keep');
    if isequal(answer,'Change'),  
        oldrad=Pro.radius;Pro.radius=0;
        while ~isempty(Pro.radius)&&~(Pro.radius>0),
            radstr=inputdlg('Specify tree radius (in m)','Input radius');
            Pro.radius=str2num(radstr{1});
            if Pro.radius>0, 
                N.elec=N.elec*Pro.radius/oldrad; 
                N.k=kfactorana(N,ana);
            end
        end
    end
    if isequal(answer,'Load File'),
        deffile=[pname filesep '*.rad'];
        [fn,pn]=uigetfile({'*.rad','radius Files(*.rad)';'*.*','All files(*.*)'},...
            'Load radius file',deffile);
        if ~isequal(fn,0), 
            Pro.radfile=fullfile(pn,fn);
            fid=fopen(Pro.radfile,'r');
            Pro.radius=fscanf(fid,'%f');
            fclose(fid);
            if max(Pro.radius)>5, Pro.radius=Pro.radius/100; end
            nel=size(N.elec,1);
            while length(Pro.radius)<nel, Pro.radius(end+1)=Pro.radius(end); end
            for i=0:nel-1, %Süden=1, dann Uhrzeigersinn
                N.elec(i+1,1)=-Pro.radius(i+1)*sin(i/nel*2*pi)*di;
                N.elec(i+1,2)=-Pro.radius(i+1)*cos(i/nel*2*pi);
            end
        end
    end
end
Pro.name=inputdlg('Type in project name:','Project Name',1,...
    {strrep(strrep(ffn,' ','_'),'+','_')});
if isempty(Pro.name), return; end
Pro.name=strrep(Pro.name{1},' ','_');
Pro.dirname=strcat(pname,Pro.name);
try,
    cd(pname);
    dos(['mkdir ' Pro.name]);
    cd(Pro.name);
    dos('mkdir tmp');
    appendmessage(handles,sprintf('Created project %s in %s',Pro.name,Pro.dirname));
    cd(progdir);
    if isequal(ffe,'.ddz')||isequal(ffe,'.tx0')||isequal(ffe,'.dat'),
        Pro.datfile=[Pro.name '.ohm'];
    end
    saveinv2dfile([Pro.dirname filesep Pro.datfile],N,1,'# x y');
    cd(Pro.dirname);
    dos(['copy /y "' fullfile(pname,fname) '" "' Pro.dirname '"']);
    cd(progdir);
    Pro.profile=[Pro.name '.pro'];
    set(handles.figure1,'Name',['DC2dTree - Project ' Pro.name]);
catch,
    disp('Catched in new profile');
    disp(lasterr);
    cd(progdir);
end
Pro.profile=strrep(Pro.profile,[Pro.dirname filesep],'');
dc2dtree('savepro_Callback',gcbo,[],guidata(gcbo));
set(handles.prepro,'Enable','On');
%if erkenne(kreis)
rads=sqrt(sum((N.elec-repmat(mean(N.elec),size(N.elec,1),1)).^2,2));
if exist(fullfile('tmp','sensMat.mat'),'file')&&exist(fullfile('tmp','meshPara.bms'),'file')&&...
        (size(N.elec,1)==24)&&(length(unique(round(rads*500)))==1),
    ff=msgbox('Circle form detected! Taking 24-electrode files');
    pause(1.0);
    if ishandle(ff), close(ff); end
%     P=loadprimpot;%im Moment nicht so richtig gebraucht
    appendmessage(handles,'Standard Circle Prepare');
    Mesh=loadmesh('tmp\meshPara.bms')
    MeshP=loadmesh('tmp\meshPrim.bms');
    MeshS=loadmesh('tmp\meshSec.bms');
    Spp=loadsens;Npp=readinv2dfile('24.data',1);
%     Npp.k=kfactorana(Npp,ana);
%     am=max(Npp.a,Npp.m)-min(Npp.a,Npp.m);
%     anan=ana(4:4:end);anan=[anan;anan(end-1:-1:1)];
%     Npp.k=1./anan(am);
    for i=1:size(Spp,1), Spp(i,:)=Spp(i,:)/Npp.k(i); end
    S=zeros(length(N.a),size(Spp,2));
    am=[Npp.a(:) Npp.m(:)];
    for i=1:length(N.a),
      S(i,:)=Spp(find(ismember(am,sort([N.a(i) N.m(i)]),'rows')),:)-...
          Spp(find(ismember(am,sort([N.a(i) N.n(i)]),'rows')),:)-...
          Spp(find(ismember(am,sort([N.b(i) N.m(i)]),'rows')),:)+...
          Spp(find(ismember(am,sort([N.b(i) N.n(i)]),'rows')),:);
      S(i,:)=S(i,:)*N.k(i);
    end
%     Pro.radius=round(mean(rads)*500)/500;
    Pro.radius=mean(rads);
%     save orient Mesh N Npp S
    [mid,ang,ori]=getdatpar(N); %mid should be [0 0]
    [midP,angP,oriP]=getdatpar(Npp);
    phi=(ang-angP)/size(N.elec,1)*2*pi;
    twister=[cos(phi) -ori*oriP*sin(phi);ori*oriP*sin(phi) cos(phi)];
    Mesh.node=Mesh.node*Pro.radius*twister;
    MeshS.node=MeshS.node*Pro.radius*twister;
    MeshP.node=MeshP.node*Pro.radius*twister;
    cd(Pro.dirname);
    savemesh(Mesh,'tmp\meshPara.bms');
    savemesh(MeshP,'tmp\meshPrim.bms');
    savemesh(MeshS,'tmp\meshSec.bms');
    savesensmat(S,'tmp\sensMat.mat');
    cd(progdir);
    [fpath,fname,fext]=fileparts(Pro.datfile);
    Pro.datafile=[fname '.data'];
    D=N;if isfield(D,'rho'), D=rmfield(D,'rho'); end
    if isfield(D,'topo'), D=rmfield(D,'topo'); end
    D.konf=D.k;
    saveinv2dfile(fullfile(Pro.dirname,Pro.datafile),D);
    set(handles.run,'Enable','On');
end
% end

% --------------------------------------------------------------------
function savepro_Callback(hObject, eventdata, handles)
global Pro progdir
if ~isfield(Pro,'dirname'), return; end
cd(Pro.dirname);
if ~isfield(Pro,'profile'), return; end
if get(handles.robust,'Value'), Pro.robust=1; end
if get(handles.blockymodel,'Value'), Pro.blockymodel=1; end
fid=fopen(Pro.profile,'w');
fprintf(fid,'DC2dTree Project File\n');
fprintf(fid,'---------------------\n');
ff=fieldnames(Pro);
for i=1:length(ff),
    nn=getfield(Pro,ff{i});
    if ~isnumeric(nn)||length(nn)<25, %
        fprintf(fid,upper(ff{i}));fprintf(fid,'=');
        if ischar(nn), fprintf(fid,'%s',nn); end
        if isnumeric(nn), 
            aa=num2strcell(nn);
            fprintf(fid,'%s ',aa{:}); 
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);
cd(progdir);
% end

% --------------------------------------------------------------------
function loadpro_Callback(hObject, eventdata, handles)
global N Pro Mesh progdir malstat ana S
infile='*.pro';
if isfield(Pro,'dirname'), infile=[Pro.dirname filesep '..' filesep '*.pro']; end
[fname,pname]=uigetfile({'*.pro','Project Files(*.pro)';...
        '*.*','All files(*.*)'},'Load data file',infile);
if isequal(fname,0), return; end
set(handles.prepro,'Enable','Off');
cd(pname);
Pro=[];Pro.radius=1;
Pro.profile=fullfile(pname,fname);
[pn,Pro.name,en]=fileparts(Pro.profile);
Pro.dirname=fileparts(Pro.profile);
set(handles.message,'String',{'DC2dTree - Tree Tomography';'-----------------------------------'});
appendmessage(handles,['Loading project file ' Pro.profile]);
fid=fopen(Pro.profile);
zeile=fgetl(fid);
while ischar(zeile),
    fi=findstr(zeile,'=');
    if fi,
        SW=zeile(1:fi-1);%fprintf('%s\n',SW);
        AR=zeile(fi+1:end);
        switch SW,
            case 'DATFILE',
                Pro.datfile=AR;
                if ~exist(Pro.datfile,'file'), 
                    uiwait(errordlg(['Dat file ' Pro.datfile ' not found!'],'File not found'));
                    fclose(fid);return;
                end
                appendmessage(handles,['Datfile = ' Pro.datfile]);
                N=readinv2dfile(Pro.datfile,1);
                if isfield(N,'rho'), Pro.rho=N.rho; end
            case 'IPFILE',
                Pro.ipfile=AR;
            case 'DATAFILE',
                Pro.datafile=AR;
                if ~exist(Pro.datafile,'file'), 
                    uiwait(errordlg(['Data file ' Pro.datafile ' not found!'],'File not found'));
                    fclose(fid);return;
                end
                appendmessage(handles,['Datafile = ' Pro.datafile]);
            case 'DDZFILE',
                Pro.ddzfile=AR;
                appendmessage(handles,['original Datfile = ' Pro.ddzfile]);
            case 'RESFILE',
                Pro.resfile=AR;
                if ~exist(Pro.resfile,'file'), 
                    uiwait(errordlg(['Result file ' Pro.resfile ' not found!'],'File not found'));
                    fclose(fid);return;
                end
                appendmessage(handles,['Result file = ' Pro.resfile]);
            case 'PHASEFILE',
                Pro.phasefile=AR;
                if ~exist(Pro.phasefile,'file'), 
                    uiwait(errordlg(['Phase model file ' Pro.phasefile ' not found!'],'File not found'));
                    fclose(fid);return;
                end
                appendmessage(handles,['Phase model file = ' Pro.phasefile]);
            case 'MORFILE',
                Pro.morfile=AR;
                if ~exist(Pro.morfile,'file'), 
                    uiwait(errordlg(['Model response file ' Pro.morfile ' not found!'],'File not found'));
                    fclose(fid);return;
                end
                appendmessage(handles,['Model response file = ' Pro.morfile]);
            case 'RADIUS',
                Pro.radius=str2num(AR);
                if length(Pro.radius)==1,
                    msg=sprintf('Radius = %.2fm',Pro.radius);
                else
                    msg=sprintf('Radius min=%.2fm,max=%.2fm',min(Pro.radius),max(Pro.radius));
                end
                appendmessage(handles,msg);
            case 'LAMBDA',
                Pro.lambda=str2num(AR);
                if ~isnumeric(Pro.lambda), Pro.lambda=30; end
                if Pro.lambda>1000, Pro.lambda=1000; end
                if Pro.lambda<1, Pro.lambda=1; end
                appendmessage(handles,['Regularization strength = ' num2str(Pro.lambda,'%.1f')]);
                set(handles.lamtext,'String',num2str(Pro.lambda,'%.1f'));
                dc2dtree('lamtext_Callback',gcbo,[],guidata(gcbo))
                cd(pname);
            case 'ROBUST',
                Pro.robust=str2num(AR);
                set(handles.robust,'Value',(Pro.robust>0));
            case 'BLOCKYMODEL',
                Pro.blockymodel=str2num(AR);
                set(handles.blockymodel,'Value',(Pro.blockymodel>0));
        end
    end
    zeile=fgetl(fid);
end
fclose(fid);
if exist(Pro.datfile,'file'),
  N=readinv2dfile(Pro.datfile,1);
  appendmessage(handles,sprintf('Reading dat file %s',Pro.datfile));
  appendmessage(handles,sprintf('%d data with %d electrodes',length(N.a),size(N.elec,1)));
  N.k=kfactorana(N,ana);
  if ~isfield(N,'rho'), N.rho=N.r; end
  N.r=abs(N.rho.*N.k);
end
if exist('tmp\sensMat.mat','file'), S=loadsens; end
if 0&&isfield(Pro,'datafile')&&exist(Pro.datafile,'file'),
    NN=readinv2dfile(Pro.datafile,1);
    appendmessage(handles,sprintf('Reading data file %s',Pro.datafile));
    if ~isfield(N,'rho'), N.rho=N.r; end
    if isfield(NN,'ip'), N.ip=NN.ip; end
%     N.r=NN.rho;N.k=NN.err*100;
end
if isfield(Pro,'ipfile')&&exist(Pro.ipfile,'file'),
    NN=readinv2dfile(Pro.ipfile,1);
    if isfield(NN,'ip'), N.ip=NN.ip; end
end
% if isfield(Pro,'ddzfile')&&exist(Pro.ddzfile,'file'),
%     NN=loadddz(Pro.ddzfile);
%     if isfield(NN,'ip'), N.ip=NN.ip; end
% end
% Mesh=loadmesh('tmp\meshPara.bms');
Mesh=[];
if exist('tmp\meshPara.bms','file'), Mesh=loadmesh('tmp\meshPara.bms'); else
    fprintf('Parameter Mesh does not exist'); end
if isfield(Pro,'resfile')&&exist(Pro.resfile), 
    fid=fopen(Pro.resfile);Pro.res=fscanf(fid,'%f');fclose(fid);
    appendmessage(handles,sprintf('Loaded model file %s',Pro.resfile));
end
if isfield(Pro,'phasefile')&&exist(Pro.phasefile), 
    fid=fopen(Pro.phasefile);Pro.phase=fscanf(fid,'%f');fclose(fid);
    appendmessage(handles,sprintf('Loaded phase file %s',Pro.phasefile));
end
if isfield(Pro,'morfile')&&exist(Pro.morfile), 
    fid=fopen(Pro.morfile);Pro.mor=fscanf(fid,'%f');fclose(fid);
    appendmessage(handles,sprintf('Loaded modelresponse file %s',Pro.morfile));
end
if exist('tmp\sensMat.mat','file')&~isempty(Mesh), 
    set(handles.run,'Enable','On'); 
    set(handles.preview,'Enable','On'); 
end
cd(progdir);
if isfield(Pro,'datafile'),
    fdf=fullfile(Pro.dirname,Pro.datafile);
    if exist(fdf,'file'), N=readinv2dfile(fdf,1); end
end
if ~isfield(N,'err')||(length(N.err)<length(N.a))||(max(N.err)<0.01), N.err=ones(size(N.a))*0.03; end
set(handles.figure1,'Name',['DC2dTree - Project ' Pro.name]);
set(handles.prepro,'Enable','On');
% if isfield(N,'r'), showpseudo(N,N.r);malstat=2; end
if ~isempty(Mesh)&&isfield(Pro,'res'), 
    dc2dtree('showmod_Callback',gcbo,[],guidata(gcbo)); 
end
% end

% --------------------------------------------------------------------
function loaddata_Callback(hObject, eventdata, handles)
global N Pro malstat ana
set(handles.message,'String',{'DC2dTree - Tree Tomography';'version 0.9.2';'------------RESISTIVITY.NET------------'});
Pro=[];
[fname,pname]=uigetfile({'*.ddz;*.ohm;*.dat;*.tx0','Known files';
    '*.ddz','Geotom Files(*.ddz)';'*.ohm','OHM Files(*.ohm)';...
    '*.dat','DAT Files(*.dat)';'*.tx0','Lippmann Files(*.tx0)';...
    '*.*','All files(*.*)'},'Load data file');
if isequal(fname,0), return; end
set(handles.prepro,'Enable','On');
set(handles.run,'Enable','Off');
set(handles.preview,'Enable','Off');
Pro=[];
[ffp,ffn,ffe]=fileparts(fname);
switch ffe,
    case '.dat',
        Pro.datfile=fullfile(pname,fname);
        N=readinv2dfile(Pro.datfile,1);
        N.k=kfactorana(N,ana);
        if ~isfield(N,'rho'), N.rho=N.r; end
        N.r=abs(N.rho.*N.k);
    case '.ohm',
        Pro.datfile=fullfile(pname,fname);
        N=readinv2dfile(Pro.datfile,1);
        N.k=kfactorana(N,ana);
        if ~isfield(N,'rho'), N.rho=N.r; end
        N.r=abs(N.rho.*N.k);
    case '.tx0',
        Pro.datfile=fullfile(pname,fname);
        N=readtx0file(Pro.datfile);
        N.k=kfactorana(N,ana);
        if ~isfield(N,'rho'), N.rho=N.r; end
        N.r=abs(N.rho.*N.k);        
    otherwise
        return;
end
appendmessage(handles,sprintf('Loaded %d data with %d electrodes',length(N.a),size(N.elec,1)));
set(handles.figure1,'Name',['DC2dTree - File ' fname]);
axes(handles.plot);
showpseudo(N,N.r);
malstat=2;
% end

% --------------------------------------------------------------------
function show_Callback(hObject, eventdata, handles)
global N Pro Mesh
onoff='off';
if isfield(N,'ip')&&any(N.ip), onoff='on'; end
set(handles.phasedata,'Enable',onoff);
onoff='off';
if isfield(Pro,'mor')&&isfield(N,'r')&&(length(Pro.mor)==length(N.r)), onoff='on'; end
set(handles.showresponse,'Enable',onoff);
set(handles.comparedata,'Enable',onoff);
onoff='off';if isfield(N,'k')&(~isempty(Mesh)), onoff='on'; end
set(handles.topoeffect,'Enable',onoff);
onoff='Off';if isfield(Pro,'phase')&&(length(Pro.phase)==length(Pro.res)), onoff='On'; end
set(handles.showphasemodel,'Enable',onoff);
% end

% --------------------------------------------------------------------
function showdata_Callback(hObject, eventdata, handles)
global N malstat
axes(handles.plot);
showpseudo(N,N.r);
malstat=2;
onoff='off';
set(handles.cmin,'Enable',onoff);
set(handles.cmintext,'Enable',onoff);
set(handles.cmax,'Enable',onoff);
set(handles.cmaxtext,'Enable',onoff);
set(handles.cauto,'Enable',onoff);
% end

% --------------------------------------------------------------------
function showresponse_Callback(hObject, eventdata, handles)
global N Pro malstat
axes(handles.plot);
showpseudo(N,Pro.mor);
malstat=3;
onoff='off';
set(handles.cmin,'Enable',onoff);
set(handles.cmintext,'Enable',onoff);
set(handles.cmax,'Enable',onoff);
set(handles.cmaxtext,'Enable',onoff);
set(handles.cauto,'Enable',onoff);
% end

% --------------------------------------------------------------------
function showmod_Callback(hObject, eventdata, handles)
global Pro Mesh N libmmfile malstat
axes(handles.plot);
if ~isfield(Pro,'res'), 
    uiwait(errordlg('No model vector loaded!'));
    return; 
end
if (length(Pro.res)==0)|(length(Pro.res)~=Mesh.ncells),
    uiwait(errordlg('Improper size of model vector!'));
    return;
end
[cmin,cmax]=tripatchmod(Mesh,Pro.res,N);
malstat=0;
global libmmfile
if ~isequal(libmmfile,4),
    xl=xlim;zl=ylim;
    set(line(xl,zl),'Color','black');
    set(line(xl,fliplr(zl)),'Color','black');
    tv=[145 144 150 140 141 154 169 223 139 140 154 171];
    tt=text(mean(xl),mean(zl),char(255-fliplr(tv)));
    set(tt,'FontSize',24,'HorizontalAlignment','center','VerticalAlignment','middle');
end
set(handles.cauto,'Enable','On','Value',1);
onoff='off';
set(handles.cmin,'Enable',onoff);
set(handles.cmintext,'Enable',onoff);
set(handles.cmax,'Enable',onoff);
set(handles.cmaxtext,'Enable',onoff);
set(handles.cmintext,'String',num2str(round(cmin*10)/10));
set(handles.cmaxtext,'String',num2str(round(cmax*10)/10));
rmed=median(Pro.res);rmin=min(Pro.res);rmax=max(Pro.res);
if rmin==rmax, return; end
if (rmin>0)&(cmin>0),
    vv1=(log10(cmin)-log10(rmin))/(log10(rmed)-log10(rmin));
    vv2=(log10(cmax)-log10(rmed))/(log10(rmax)-log10(rmed));
else
    vv1=(cmin-rmin)/(rmed-rmin);
    vv2=(cmax-rmed)/(rmax-rmed);
end
if (vv1>0)&(vv1<1), set(handles.cmin,'Value',vv1); end
if (vv2>0)&(vv2<1), set(handles.cmax,'Value',vv2); end
% end

% --- Executes during object creation, after setting all properties.
function lamslide_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% end

% --- Executes on slider movement.
function lamslide_Callback(hObject, eventdata, handles)
val=get(handles.lamslide,'Value');
lam=10^(3*(1-val));
set(handles.lamtext,'String',num2str(lam,'%.1f'));
dc2dtree('lamtext_Callback',gcbo,[],guidata(gcbo))
% end

% --- Executes during object creation, after setting all properties.
function lamtext_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% end

function lamtext_Callback(hObject, eventdata, handles)
lam=str2num(get(handles.lamtext,'String'));
set(handles.lamslide,'Val',1-log10(lam)/3);
dc2dtree('preview_Callback',gcbo,[],guidata(gcbo))
% end

% --- Executes on button press in lamdefault.
function lamdefault_Callback(hObject, eventdata, handles)
lam=30;
set(handles.lamtext,'String',num2str(lam));
set(handles.lamslide,'Val',1-log10(lam)/3);
dc2dtree('preview_Callback',gcbo,[],guidata(gcbo))
% end

% --- Executes on button press in newpro1.
function newpro1_Callback(hObject, eventdata, handles)
% end

% --------------------------------------------------------------------
function exitprogram_Callback(hObject, eventdata, handles)
for i=1:10,
    if ishandle(i), delete(i); end
end
delete(gcbf);
% end

% --------------------------------------------------------------------
function export_Callback(hObject, eventdata, handles)
global Pro N
onoff='off';
if isfield(N,'ip')&&any(N.ip), onoff='on'; end
set(handles.phasedata,'Enable',onoff);
ismod=sonoff(isfield(Pro,'res'));
set(handles.exportmod,'Enable',ismod);
set(handles.exportvtk,'Enable',ismod);
set(handles.exportphase,'Enable',sonoff(isfield(Pro,'phase')));
set(handles.probeline,'Enable',ismod);%new feature
% end

% --------------------------------------------------------------------
function exportmod_Callback(hObject, eventdata, handles)
global Pro Mesh N progdir
f=figure;
set(f,'Units','Characters','MenuBar','none','NumberTitle','off','Name','Model');
iconify(f);
po=get(f,'Position');
po1=get(handles.plot,'Position');
po(3)=po1(3)*1.3;
set(f,'Position',po);
if get(handles.cauto,'Value')==1,
    tripatchmod(Mesh,Pro.res,N);
else
    cmin=str2num(get(handles.cmintext,'String'));
    cmax=str2num(get(handles.cmaxtext,'String'));
    tripatchmod(Mesh,Pro.res,N,cmin,cmax);
end
global libmmfile
if ~isequal(libmmfile,4),
    xl=xlim;zl=ylim;
    set(line(xl,zl),'Color','black');
    set(line(xl,fliplr(zl)),'Color','black');
    tv=[145 144 150 140 141 154 169 223 139 140 154 171];
    tt=text(mean(xl),mean(zl),char(255-fliplr(tv)));
    set(tt,'FontSize',24,'HorizontalAlignment','center','VerticalAlignment','middle');
end
outfile=strcat(Pro.dirname,filesep,Pro.name,'.pdf');
[fname,pname]=uiputfile('*.pdf','Save figure as',outfile);
if ~isequal(fname,0),
    ff=strrep(fullfile(pname,fname),'.pdf','');
    epsprint(f,ff,1);
    dos([progdir filesep 'epstopdf "' ff '.eps"']);
%     delete([ff '.eps']);
    appendmessage(handles,['Exported model as pdf file ' ff '.pdf']);
end
close(f);
% end

% --------------------------------------------------------------------
function exportvtk_Callback(hObject, eventdata, handles)
global Pro Mesh libmmfile
if ~isequal(libmmfile,4), warndlg('Not available in Test version'); return; end
if ~isfield(Pro,'res'), return; end
outfile=strcat(Pro.dirname,filesep,Pro.name,'.vtk');
[fname,pname]=uiputfile('*.vtk','Save model as',outfile);
if isequal(fname,0), return; end
Mesh.cellattr=Pro.res;
if isfield(Pro,'phase'),
    savevtkmesh(Mesh,fullfile(pname,fname),Pro.phase);
else
    savevtkmesh(Mesh,fullfile(pname,fname));
end
appendmessage(handles,['Exported model as vtk file ' fullfile(pname,fname)]);
return
fid=fopen(fullfile(pname,fname),'w');
fprintf(fid,'# vtk DataFile Version 3.0\ncreated by DC2dTree\nASCII\nDATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %d double\n',Mesh.nnodes);
fprintf(fid,'%.3f\t%.3f\t%.3f\n',Mesh.node');
nel=Mesh.ncells;
fprintf(fid,'CELLS %d %d\n',nel,nel*4);
fprintf(fid,'3\t%d\t%d\t%d\n',Mesh.cell'-1);
fprintf(fid,'CELL_TYPES %d\n',nel);
fprintf(fid,'%d ',ones(nel,1)*5);
fprintf(fid,'\nCELL_DATA %d\nSCALARS log10(Resistivity) double 1\nLOOKUP_TABLE default\n',nel);
fprintf(fid,'%.3f ',log10(Pro.res));
fprintf(fid,'\nSCALARS Resistivity double 1\nLOOKUP_TABLE default\n');
fprintf(fid,'%.3f ',Pro.res);
fprintf(fid,'\nSCALARS Material double 1\nLOOKUP_TABLE default\n');
fprintf(fid,'%d ',2:nel+1);
fprintf(fid,'\n');
fclose(fid);
appendmessage(handles,['Exported model as vtk file ' fullfile(pname,fname)]);
% end

% --------------------------------------------------------------------
function exportdata_Callback(hObject, eventdata, handles)
global Pro N progdir
f=figure;
set(f,'Units','Characters','MenuBar','none','NumberTitle','off','Name','Model');
iconify(f);
po=get(f,'Position');
po1=get(handles.plot,'Position');
po(3)=po1(3)*1.3;
set(f,'Position',po);
showpseudo(N,N.r);
outfile=strcat(Pro.dirname,filesep,Pro.name,'-data.pdf');
[fname,pname]=uiputfile('*.pdf','Save figure as',outfile);
if ~isequal(fname,0),
    outfile=fullfile(pname,fname);
    appendmessage(handles,sprintf('Exported data as pdf %s',outfile));
    ff=strrep(outfile,'.pdf','');
    epsprint(f,ff,1);
    dos([progdir filesep 'epstopdf "' ff '.eps"']);
end
close(f);
% end

% --- Executes on button press in optlambda.
function optlambda_Callback(hObject, eventdata, handles)
onoff='On';
if get(handles.optlambda,'Value'), onoff='Off'; end
set(handles.lamslide,'Enable',onoff);
set(handles.lamtext,'Enable',onoff);
set(handles.preview,'Enable',onoff);
% end

% --- Executes on key press over figure1 with no controls selected.
function figure1_KeyPressFcn(hObject, eventdata, handles)
global S Mesh N Pro progdir
aa=get(gcf,'CurrentCharacter');
set(gcf,'Pointer','arrow');
switch aa,
    case '1',
        dc2dtree('onetv_Callback',gcbo,[],guidata(gcbo));
    case 'D',
        dc2dtree('showdata_Callback',gcbo,[],guidata(gcbo));
    case 'I',
        dc2dtree('phaseinv_Callback',gcbo,[],guidata(gcbo));
    case 'G',
        dc2dtree('optmesh_Callback',gcbo,[],guidata(gcbo));
    case 'L',
        dc2dtree('loadpro_Callback',gcbo,[],guidata(gcbo));
    case 'M',
        dc2dtree('showmod_Callback',gcbo,[],guidata(gcbo));
    case 'N',
        dc2dtree('newpro_Callback',gcbo,[],guidata(gcbo));
    case 'P',
        dc2dtree('prepro_Callback',gcbo,[],guidata(gcbo));
    case 'R',
        dc2dtree('run_Callback',gcbo,[],guidata(gcbo));        
    case 'S',
        if isempty(S),
            cd(Pro.dirname);S=loadsens;Mesh=loadmesh('tmp\meshPara');cd(progdir);
        end
        l=floor(rand*size(S,1))+1;axes(handles.plot);
        tripatchmod(Mesh,S(l,:),N);
        abmn=[N.a(l) N.b(l) N.m(l) N.n(l)];
        hold on;plot(N.elec(abmn,1),N.elec(abmn,2),'o');hold off
    case 'T',
        dc2dtree('topoeffect_Callback',gcbo,[],guidata(gcbo));        
    case 'X',
        dc2dtree('exitprogram_Callback',gcbo,[],guidata(gcbo));        
    case 'Z',
        dc2dtree('timeall_Callback',gcbo,[],guidata(gcbo));        
end
% end

% --------------------------------------------------------------------
function phasedata_Callback(hObject, eventdata, handles)
global N malstat
axes(handles.plot);
if ~isfield(N,'ip'), return; end
ip=N.ip;ip(find(N.ip>=30))=NaN;
showpseudo(N,ip);
malstat=4;
onoff='off';
set(handles.cmin,'Enable',onoff);
set(handles.cmintext,'Enable',onoff);
set(handles.cmax,'Enable',onoff);
set(handles.cmaxtext,'Enable',onoff);
set(handles.cauto,'Enable',onoff);
% end

% --------------------------------------------------------------------
function do_Callback(hObject, eventdata, handles)
global N Pro S Mesh
onoff='off';
if isfield(N,'r')&&(length(N.r)>0), 
    onoff='on'; end
set(handles.prepro,'Enable',sonoff);
onoff='off';
if isfield(N,'r')&&(length(N.r)==size(S,1)), 
    onoff='on'; end
set(handles.run,'Enable',sonoff);
onoff='off';
if isfield(N,'ip')&&(length(N.ip)==length(N.r))&&isfield(Pro,'res'), 
    onoff='on'; 
end
set(handles.phaseinv,'Enable',onoff);
% set(handles.batchinversion,'Enable','off'); %new feature
% end

% --------------------------------------------------------------------
function phaseinv_Callback(hObject, eventdata, handles)
global Pro N Mesh progdir malstat
if isfield(Pro,'dirname'), cd(Pro.dirname); else return; end
Pro.lambda=str2num(get(handles.lamtext,'String'));
S=loadsens;
if get(handles.blockymodel,'Value')&&exist('CRobust.mat','file'),
    fid=fopen('CRobust.mat','r');
else
    fid=fopen(['tmp' filesep 'smoothness.matrix'],'r');
end
if fid<1, uiwait(errordlg('Did not find smoothness matrix!'));return; end
mm=fscanf(fid,'%f',[3 Inf])';
fclose(fid);
cd(progdir);
C=sparse(mm(:,1)+1,mm(:,2)+1,mm(:,3))+speye(Mesh.ncells)*0.03;
D=spdiags(1./log(N.err+1),0,length(N.err),length(N.err));
set(gcf,'Pointer','watch');
if 0,
    Pro.phase=cglscdp(S,N.ip,Pro.lambda,C,D);
else
    fi=find((N.ip<0)&(abs(N.ip)<10));
    D=spdiags(1./log(N.err(fi)+1),0,length(fi),length(fi));
    if 0,
        dip=N.ip(fi);
        Pro.phase=cglscdp(S(fi,:),dip,Pro.lambda,C,D);
    else
        if strcmp(lower(Pro.ddzfile(end-3:end)),'.tx0'),
            Org=readtx0file(Pro.ddzfile);
        else
            Org=[];
        end
        dip=log(-N.ip(fi));
        Pro.phase=-exp(cglscdp(S(fi,:),dip,Pro.lambda,C,D));
        if isfield(Org,'allip')&&(size(Org.allip,1)==length(N.ip)),
            nfr=size(Org.allip,2);
            Pro.phases=repmat(Pro.phase,1,nfr);
            for i=2:nfr,
                fi=find((Org.allip(:,i)<0)&(abs(N.ip)<10));
                D=spdiags(1./log(N.err(fi)+1),0,length(fi),length(fi));
                dip=log(-Org.allip(fi,i));
                Pro.phases(:,i)=-exp(cglscdp(S(fi,:),dip,Pro.lambda,C,D));                
            end
        end
    end
end
set(gcf,'Pointer','arrow');
Pro.phasefile=[Pro.name '.phase'];
cd(Pro.dirname);
fid=fopen(Pro.phasefile,'w');fprintf(fid,'%g\n',Pro.phase(:));fclose(fid);
cd(progdir);
malstat=1;
tripatchmod(Mesh,Pro.phase,N);
onoff='off';
set(handles.cmin,'Enable',onoff);
set(handles.cmintext,'Enable',onoff);
set(handles.cmax,'Enable',onoff);
set(handles.cmaxtext,'Enable',onoff);
set(handles.cauto,'Enable',onoff);
dc2dtree('savepro_Callback',gcbo,[],guidata(gcbo));
% end

% --- Executes during object creation, after setting all properties.
function cmin_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% end

% --- Executes on slider movement.
function cmintext_Callback(hObject, eventdata, handles)
global Pro malstat
if (malstat==1)&&isfield(Pro,'phase'), field=Pro.phase;
else field=Pro.res; end
rmed=median(field);rmin=min(field);rmax=max(field);
cmin=str2num(get(handles.cmintext,'String'));
if rmin>0,
    vv=(log10(cmin)-log10(rmin))/(log10(rmed)-log10(rmin));
else
    vv=(cmin-rmin)/(rmed-rmin);
end
if (vv>0)&(vv<1), set(handles.cmin,'Value',vv); end
drawnow;
dc2dtree('cdraw_Callback',gcbo,[],guidata(gcbo))
% end

% --- Executes on slider movement.
function cmin_Callback(hObject, eventdata, handles)
global Pro malstat
if (malstat==1)&&isfield(Pro,'phase'), field=Pro.phase;
else field=Pro.res; end
rmed=median(field);rmin=min(field);rmax=max(field);
vv=get(handles.cmin,'Value');
if rmin>0,
    cmin=10^(vv*log10(rmed/rmin))*rmin;
else
    cmin=vv*(rmed-rmin)+rmin;
end
set(handles.cmintext,'String',num2str(round(cmin*10)/10));
dc2dtree('cdraw_Callback',gcbo,[],guidata(gcbo))
% end

% --- Executes during object creation, after setting all properties.
function cmax_CreateFcn(hObject, eventdata, handles)
% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% end

% --- Executes on slider movement.
function cmaxtext_Callback(hObject, eventdata, handles)
global Pro malstat
if (malstat==1)&&isfield(Pro,'phase'), field=Pro.phase;
else field=Pro.res; end
rmed=median(field);rmax=max(field);rmin=min(field);
cmax=str2num(get(handles.cmaxtext,'String'));
if rmin>0,
    vv=(log10(cmax)-log10(rmed))/(log10(rmax)-log10(rmed));
else
    vv=(cmax-rmed)/(rmax-rmed);
end
if (vv>0)&(vv<1), set(handles.cmax,'Value',vv); end
drawnow;
dc2dtree('cdraw_Callback',gcbo,[],guidata(gcbo))
% end

% --- Executes on slider movement.
function cmax_Callback(hObject, eventdata, handles)
global Pro malstat
if (malstat==1)&&isfield(Pro,'phase'), field=Pro.phase;
else field=Pro.res; end
rmed=median(field);rmin=min(field);rmax=max(field);
vv=get(handles.cmax,'Value');
if rmin>0,
    cmax=10^(vv*log10(rmax/rmed))*rmed;
else
    cmax=vv*(rmax-rmed)+rmed;
end
set(handles.cmaxtext,'String',num2str(round(cmax*10)/10));
dc2dtree('cdraw_Callback',gcbo,[],guidata(gcbo))
% end

% --------------------------------------------------------------------
function optmesh_Callback(hObject, eventdata, handles)
hh=meshopts;
iconify(hh);
uiwait(hh);
% end

% --------------------------------------------------------------------
function optinv_Callback(hObject, eventdata, handles)
global Pro Invopts
hh=invopts;
iconify(hh);
uiwait(hh);
set(handles.blockymodel,'Value',fvalue(Invopts,'blockymodel'));
set(handles.robust,'Value',fvalue(Invopts,'robustdata'));
set(handles.optlambda,'Value',fvalue(Invopts,'optlambda'));
% end

% --------------------------------------------------------------------
function options_Callback(hObject, eventdata, handles)
% end

% --------------------------------------------------------------------
function Help_Callback(hObject, eventdata, handles)
global progdir Mesh Pro
cd(progdir);
% disp(Pro);
% disp(Mesh);
% end

% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
% hh=about;
hh=msgbox({'DC2dTree';'Impedance Tomography';'on trees with tree shapes';...
        'version 0.9.2';'';'Authos: C. Rücker & T. Günther';'resistivity.net productions'},'About');
iconify(hh);
uiwait(hh);
% end

% --- Executes on button press in cauto.
function cauto_Callback(hObject, eventdata, handles)
global Pro malstat
onoff='on';
if get(handles.cauto,'Value')==1, 
    onoff='off';
    switch malstat,
        case 1,
            dc2dtree('showphasemodel_Callback',gcbo,[],guidata(gcbo))
        otherwise
            dc2dtree('showmod_Callback',gcbo,[],guidata(gcbo))
    end
end
set(handles.cmin,'Enable',onoff);
set(handles.cmintext,'Enable',onoff);
set(handles.cmax,'Enable',onoff);
set(handles.cmaxtext,'Enable',onoff);
% end

function cdraw_Callback(hObject, eventdata, handles)
global Pro Mesh N libmmfile malstat
axes(handles.plot);
if ~isfield(Pro,'res'), 
    uiwait(errordlg('No model vector loaded!'));
    return; 
end
if (isempty(Pro.res))||(length(Pro.res)~=Mesh.ncells),
    uiwait(errordlg('Improper size of model vector!'));
    return;
end
cmin=str2num(get(handles.cmintext,'String'));
cmax=str2num(get(handles.cmaxtext,'String'));
switch malstat,
    case 1, 
        tripatchmod(Mesh,Pro.phase,N,cmin,cmax);
        if isfield(Pro,'freq')&&(length(Pro.freq)>=1), 
            xl=xlim;yl=ylim;t=text(xl(1),yl(2),sprintf(' f=%gHz',Pro.freq(1))); 
            set(t,'VerticalAlignment','top');
        end
        if isfield(Pro,'phases')&&(size(Pro.phase,1)==length(Pro.phase))&&(size(Pro.phases,2)>1), %SIP phase
            nf=size(Pro.phases,2);
            f=figure(1);iconify(f);
            set(f,'Units','Pixels','MenuBar','none','NumberTitle','off','Name','SIP phases');
            nu=fix(sqrt(nf-1));nv=fix((nf-2)/nu)+1;
            for i=2:nf, 
                subplot(nu,nv,i-1);tripatchmod(Mesh,Pro.phases(:,i),N,cmin,cmax); 
                if isfield(Pro,'freq')&&(length(Pro.freq)>=i), title(sprintf('f=%gHz',Pro.freq(i))); end
            end
        end        
otherwise
    tripatchmod(Mesh,Pro.res,N,cmin,cmax);
    global libmmfile
    if ~isequal(libmmfile,4),
        xl=xlim;zl=ylim;
        set(line(xl,zl),'Color','black');
        set(line(xl,fliplr(zl)),'Color','black');
        tv=[145 144 150 140 141 154 169 223 139 140 154 171];
        tt=text(mean(xl),mean(zl),char(255-fliplr(tv)));
        set(tt,'FontSize',24,'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end
% end

% --------------------------------------------------------------------
function topoeffect_Callback(hObject, eventdata, handles)
global Pro N ana
kfalsch=kfactorana(N,ana);
topoeff=abs(kfalsch./N.k);
f=figure(1);
set(f,'Units','Pixels','MenuBar','none','NumberTitle','off','Name','Topographic correction');
iconify(f);
po=get(f,'Position');po(3)=840;po(4)=400;
set(f,'Position',po);
rawdata=N.r.*topoeff;
cmax=max([rawdata(:);N.r(:)]);
cmin=min([rawdata(:);N.r(:)]);
subplot(1,3,1);
showpseudo(N,rawdata,cmin,cmax);
title('Raw data');
subplot(1,3,3);
showpseudo(N,N.r,cmin,cmax);
title('Corrected data');
subplot(1,3,2);
cmax=max(max(topoeff),max(1./topoeff));cmin=1/cmax;
showpseudo(N,topoeff,cmin,cmax);
title('Topography effect');
% onoff='off';
% set(handles.cmin,'Enable',onoff);
% set(handles.cmintext,'Enable',onoff);
% set(handles.cmax,'Enable',onoff);
% set(handles.cmaxtext,'Enable',onoff);
% set(handles.cauto,'Enable',onoff);
% end

% --------------------------------------------------------------------
function comparedata_Callback(hObject, eventdata, handles)
global Pro N
if isfield(Pro,'mor')&&(length(N.r)==length(Pro.mor)),
    f=figure(1);
    set(f,'Units','Pixels','MenuBar','none','NumberTitle','off','Name','Data, model response and misfit');
    iconify(f);
    po=get(f,'Position');po(3)=840;po(4)=400;
    set(f,'Position',po);
    cmin=min(N.r);cmax=max(N.r);
    subplot(1,3,1);
    showpseudo(N,N.r,cmin,cmax);
    title('Measured data');
    subplot(1,3,2);
    showpseudo(N,Pro.mor,cmin,cmax);
    title('Model response');
    misfit=(Pro.mor./N.r-1)*100;
    subplot(1,3,3);
    cmax=max(abs(misfit));cmin=-cmax;
    showpseudo(N,misfit,cmin,cmax);
    title('Misfit in %');
end
% end

% --------------------------------------------------------------------
function onetv_Callback(hObject, eventdata, handles)
global Pro progdir N Mesh
if isfield(Pro,'dirname'), cd(Pro.dirname); else return; end
Pro.lambda=str2num(get(handles.lamtext,'String'));
S=loadsens;
fid=fopen(['tmp' filesep 'smoothness.matrix'],'r');
if fid<1, uiwait(errordlg('Did not find smoothness matrix!'));return; end
mm=fscanf(fid,'%f',[3 Inf])';
fclose(fid);
cd(progdir);
ep=0.3;ep2=ep^2;ep4=ep2^2;nm=size(S,2);
C=sparse(mm(:,1)+1,mm(:,2)+1,mm(:,3));
for i=1:nm, C(i,i)=0; end
for i=1:size(mm,1),
    ii=mm(i,1)+1;jj=mm(i,2)+1;
    if ii>jj,
        dm=log(Pro.res(ii))-log(Pro.res(jj));
        h1=ep2*(ep2-3*dm^4)./(dm.^2+ep2)^4;
        C(ii,jj)=h1;C(jj,ii)=h1;
        C(ii,ii)=C(ii,ii)-h1;
        C(jj,jj)=C(jj,jj)-h1;
    end
end
f=C*log(Pro.res);lam=0.1;
dwid=(S'*S+lam*C)\(S'*(log(Pro.mor)-log(N.r))-lam*f);
% minmax(dwid)
axes(handles.plot);
% tripatchmod(Mesh,dwid,N);
newres=Pro.res.*exp(dwid);
tripatchmod(Mesh,newres,N,min(newres),max(newres));
newmor=Pro.mor.*exp(S*dwid);
[rms(N.r,Pro.mor) rms(N.r,newmor)]
Pro.mor=newmor;Pro.res=newres;
% end

% --------------------------------------------------------------------
function exportip_Callback(hObject, eventdata, handles)
global Pro N progdir
f=figure;
set(f,'Units','Characters','MenuBar','none','NumberTitle','off','Name','Model');
iconify(f);
po=get(f,'Position');
po1=get(handles.plot,'Position');
po(3)=po1(3)*1.3;
set(f,'Position',po);
ip=N.ip;ip(find(N.ip>=0))=NaN;
showpseudo(N,ip);
outfile=strcat(Pro.dirname,filesep,Pro.name,'-ipdata.pdf');
[fname,pname]=uiputfile('*.pdf','Save figure as',outfile);
if ~isequal(fname,0),
    outfile=fullfile(pname,fname);
    appendmessage(handles,sprintf('Exported data as pdf %s',outfile));
    ff=strrep(outfile,'.pdf','');
    epsprint(f,ff,1);
    dos([progdir filesep 'epstopdf "' ff '.eps"']);
end
close(f);
% end

% --------------------------------------------------------------------
function exportphase_Callback(hObject, eventdata, handles)
global Pro Mesh N progdir malstat
f=figure;
set(f,'Units','Characters','MenuBar','none','NumberTitle','off','Name','Model');
iconify(f);
po=get(f,'Position');
po1=get(handles.plot,'Position');
po(3)=po1(3)*1.3;
set(f,'Position',po);
if (get(handles.cauto,'Value')==0)&&(malstat==1),
    cmin=str2num(get(handles.cmintext,'String'));
    cmax=str2num(get(handles.cmaxtext,'String'));
    tripatchmod(Mesh,Pro.phase,N,cmin,cmax);
else
    tripatchmod(Mesh,Pro.phase,N);
end
% tripatchmod(Mesh,Pro.phase,N);
global libmmfile
if ~isequal(libmmfile,4),
    xl=xlim;zl=ylim;
    set(line(xl,zl),'Color','black');
    set(line(xl,fliplr(zl)),'Color','black');
    tv=[145 144 150 140 141 154 169 223 139 140 154 171];
    tt=text(mean(xl),mean(zl),char(255-fliplr(tv)));
    set(tt,'FontSize',24,'HorizontalAlignment','center','VerticalAlignment','middle');
end
outfile=strcat(Pro.dirname,filesep,Pro.name,'-phase.pdf');
[fname,pname]=uiputfile('*.pdf','Save figure as',outfile);
if ~isequal(fname,0),
    ff=strrep(fullfile(pname,fname),'.pdf','');
    epsprint(f,ff,1);
    dos([progdir filesep 'epstopdf "' ff '.eps"']);
%     delete([ff '.eps']);
    appendmessage(handles,['Exported phase as pdf file ' ff '.pdf']);
end
close(f);
% end

% --- Executes on button press in preview.
function preview_Callback(hObject, eventdata, handles)
global S Pro N progdir Mesh malstat
if get(handles.preview,'Value'),
    cd(Pro.dirname);
    if ~isequal(size(S,1),length(N.r)), S=loadsens; end
    fid=fopen(['tmp' filesep 'smoothness.matrix'],'r');
    if fid<1, uiwait(errordlg('Did not find smoothness matrix!'));return; end
    mm=fscanf(fid,'%f',[3 Inf])';
    fclose(fid);
    cd(progdir);
    C=sparse(mm(:,1)+1,mm(:,2)+1,mm(:,3))+speye(Mesh.ncells)*0.03;
    D=spdiags(1./log(N.err+1),0,length(N.err),length(N.err));
    rho0=median(N.r);
    dr=log(N.r)-log(rho0);
    lam=str2num(get(handles.lamtext,'String'));
    set(gcf,'Pointer','watch');
    tmpmodel=rho0.*exp(cglscdp(S,dr,lam,C,D));
    set(gcf,'Pointer','arrow');
    axes(handles.plot);
    malstat=0;
    tripatchmod(Mesh,tmpmodel,N);
end
cd(progdir);
% end

% --- Executes on button press in robust.
function robust_Callback(hObject, eventdata, handles)
% end

% --- Executes on button press in blockymodel.
function blockymodel_Callback(hObject, eventdata, handles)
% end

% --------------------------------------------------------------------
function showphasemodel_Callback(hObject, eventdata, handles)
global Pro Mesh N libmmfile malstat
axes(handles.plot);
if ~isfield(Pro,'phase'), 
    uiwait(errordlg('No phase model vector loaded!'));
    return; 
end
if (length(Pro.phase)==0)|(length(Pro.phase)~=Mesh.ncells),
    uiwait(errordlg('Improper size of model vector!'));
    return;
end
[cmin,cmax]=tripatchmod(Mesh,Pro.phase,N);
malstat=1;
global libmmfile
if ~isequal(libmmfile,4),
    xl=xlim;zl=ylim;
    set(line(xl,zl),'Color','black');
    set(line(xl,fliplr(zl)),'Color','black');
    tv=[145 144 150 140 141 154 169 223 139 140 154 171];
    tt=text(mean(xl),mean(zl),char(255-fliplr(tv)));
    set(tt,'FontSize',24,'HorizontalAlignment','center','VerticalAlignment','middle');
end
if isfield(Pro,'phases')&&(size(Pro.phase,1)==length(Pro.phase))&&(size(Pro.phases,2)>1), %SIP phase
    if isfield(Pro,'freq')&&(length(Pro.freq)>=1), 
        xl=xlim;yl=ylim;t=text(xl(1),yl(2),sprintf(' f=%gHz',Pro.freq(1))); 
        set(t,'VerticalAlignment','top');
    end
    nf=size(Pro.phases,2);
    f=figure(1);iconify(f);
    set(f,'Units','Pixels','MenuBar','none','NumberTitle','off','Name','SIP phases');
    nu=fix(sqrt(nf-1));nv=fix((nf-2)/nu)+1;
    for i=2:nf, 
        subplot(nu,nv,i-1);tripatchmod(Mesh,Pro.phases(:,i)); 
        if isfield(Pro,'freq')&&(length(Pro.freq)>=i), title(sprintf('f=%gHz',Pro.freq(i))); end
    end
end
set(handles.cauto,'Enable','On','Value',1);
onoff='off';
set(handles.cmin,'Enable',onoff);
set(handles.cmintext,'Enable',onoff);
set(handles.cmax,'Enable',onoff);
set(handles.cmaxtext,'Enable',onoff);
set(handles.cmintext,'String',num2str(round(cmin*10)/10));
set(handles.cmaxtext,'String',num2str(round(cmax*10)/10));
rmed=median(Pro.phase);rmin=min(Pro.phase);rmax=max(Pro.phase);
if rmin==rmax, return; end
%% if log!!!
vv=(log10(cmin)-log10(rmin))/(log10(rmed)-log10(rmin));
if (vv>0)&(vv<1), set(handles.cmin,'Value',vv); end
vv=(log10(cmax)-log10(rmed))/(log10(rmax)-log10(rmed));
if (vv>0)&(vv<1), set(handles.cmax,'Value',vv); end

function appendmessage(handles,msgtext)
if nargin<2, display('Wrong appendmessage');return; end
if isfield(handles,'message')&&ishandle(handles.message),
    mess=get(handles.message,'String');
    ma=36;
    while length(msgtext)>ma,
        while(double(msgtext)<58)&(ma>30), ma=ma-1; end
        mess{end+1}=msgtext(1:ma);
        msgtext(1:ma-2)='';
        msgtext(1:2)='  ';
        ma=36;
    end
    if length(msgtext)>0, mess{end+1}=msgtext; end
    set(handles.message,'String',mess);
else
    display(msgtext);
end
% end

% --------------------------------------------------------------------
function xyzfile_Callback(hObject, eventdata, handles)
global Pro Mesh
A=zeros(Mesh.ncells,2);ss='%.4f\t%.4f';
for i=1:Mesh.ncells, 
   A(i,:)=mean(Mesh.node(Mesh.cell(i,:),1:2));
end
if isfield(Pro,'res')&&(length(Pro.res)==Mesh.ncells),
    A(:,end+1)=Pro.res(:);ss=[ss '\t%.4f']; end
if isfield(Pro,'phase')&&(length(Pro.phase)==Mesh.ncells),
    A(:,end+1)=Pro.phase(:);ss=[ss '\t%.4f']; end
outfile=strcat(Pro.dirname,filesep,Pro.name,'.xyz');
[fname,pname]=uiputfile('*.xyz','Save  as',outfile);
if isequal(fname,0), return; end
% save(fullfile(pname,fname),'A','-ascii');
fid=fopen(fullfile(pname,fname),'w');
fprintf(fid,[ss '\r\n'],A');
fclose(fid);
% end

% --------------------------------------------------------------------
function pickvalue_Callback(hObject, eventdata, handles)
global Mesh Pro
cellmids=zeros(Mesh.ncells,2);
for i=1:Mesh.ncells, cellmids(i,:)=mean(Mesh.node(Mesh.cell(i,:),:)); end
[x,y]=ginput(1);
di=sum((cellmids-repmat([x y],Mesh.ncells,1)).^2,2);
[midi,nn]=min(di);
aa=sprintf('Cell No. %d: x=%.3fm y=%.3fm',nn,cellmids(nn,1),cellmids(nn,2));
if isfield(Pro,'res')&&(length(Pro.res)==Mesh.ncells),
    aa=sprintf('%s rho=%.1f Ohmm',aa,Pro.res(nn)); end
if isfield(Pro,'phase')&&(length(Pro.phase)==Mesh.ncells),
    aa=sprintf('%s phi=%.1f mrad',aa,Pro.phase(nn)); end
nnn=Mesh.cell(nn,:);nnn(end+1)=nnn(1);
hold on;line(Mesh.node(nnn,1),Mesh.node(nnn,2),'Color','white');hold off
if isfield(Pro,'phases')&&(size(Pro.phases,1)>=nn),
    f=figure(3);set(f,'MenuBar','none','NumberTitle','off','Name','Picked cell');
    iconify(f);
    plot(1:size(Pro.phases,2),Pro.phases(nn,:),'x-');ylabel('phi in mrad');
    title(aa);
else
    iconify(msgbox(aa,'Picked cell'));
end
% end

function erg=sonoff(value)
erg='Off';
if (nargin>0)&&isnumeric(value)&&(value>0), erg='On'; end
if (nargin>0)&&islogical(value)&&(value), erg='On'; end
% end

function erg=fonoff(Struct,field)
erg='Off';
if (nargin>1)&&isfield(Struct,field)&&(getfield(Struct,field)>0), erg='On'; end
% end

function erg=fvalue(Struct,field)
erg=0;
if (nargin>1)&&isfield(Struct,field)&&(getfield(Struct,field)>0), erg=1; end
% end

% --------------------------------------------------------------------
function timelapse_Callback(hObject, eventdata, handles)
global TL N
sonoff='Off';
if isfield(TL,'data')&&(size(TL.data,1)==length(N.r))&&(size(TL.data,2)>1), sonoff='On'; end
set(handles.timedata,'Enable',sonoff);
set(handles.timeinv,'Enable',sonoff);
sonoff='Off';
if isfield(TL,'deltamodel')&&(~isempty(TL.deltamodel)), sonoff='On'; end
set(handles.timeexport,'Enable',sonoff);
set(handles.timepickcell,'Enable',sonoff);%new feature
% end 

% --------------------------------------------------------------------
function timeload_Callback(hObject, eventdata, handles)
global N TL Pro
set(gcf,'Pointer','watch');
sf=strfind(Pro.dirname,filesep);
alltypes='*.ddz;*.ohm;*.dat';
TL=[];
[TL.first,TL.dir]=uigetfile(alltypes,'Select reference file',fullfile(Pro.dirname(1:sf(end)-1),alltypes));
if ~ischar(TL.first), return; end
[pp,ff,ee]=fileparts(TL.first);isddz=0;
switch ee,
case '.ddz',
    N0=loadddz(fullfile(TL.dir,TL.first));isddz=1;
otherwise
    N0=sort2delecs(readinv2dfile(fullfile(TL.dir,TL.first),1));
end
if ~isfield(N0,'r'),
    if ~isfield(N0,'rho0'),
        if isfield(N0,'u')&isfield(N0,'i'), rho0=N0.u./N0.i; end    
    end
    if isfield(N0,'rho')&isfield(N,'k'), N0.r=N0.rho(:).*N.k(:); end    
end
aa=dir(fullfile(TL.dir,['*' ee]));
dd=zeros(length(aa),1);
for i=1:length(aa), dd(i)=datenum(aa(i).date); end
[dds,so]=sort(dd);
if isempty(aa), uiwait(errordlg('Did not find files!'));return; end
j=0;TL.names={};TL.data=zeros(length(N.a),1);
abmn0=[N0.a N0.b N0.m N0.n];
for i=1:length(aa),
    if isddz, NN=loadddz(fullfile(TL.dir,aa(so(i)).name));
    else NN=sort2delecs(readinv2dfile(fullfile(TL.dir,aa(so(i)).name),1)); end
    if length(N0.a)==length(NN.a),
        if isfield(NN,'rho'), rho=NN.rho; else rho=NN.u./NN.i; end
        rhoa=rho.*N.k;
         if min(rhoa)>0, %max(abs(deltadata))<0.5,
            deltadata=log(rhoa)-log(N0.r);
            j=j+1;
            TL.data(:,j)=rho.*N.k;
%             data(:,j)=deltadata;
            TL.names{j}=aa(so(i)).date(end-7:end);%strrep(aa(i).name,'.ddz','');
            TL.dates(j)=datenum(aa(so(i)).date);
        else
            fprintf('Neglecting measurement %d: %s\n',i,aa(i).name);
        end
    end
end
set(gcf,'Pointer','arrow');
% end

% --------------------------------------------------------------------
function timedata_Callback(hObject, eventdata, handles)
global TL
figure(1);
set(1,'Units','Characters','MenuBar','none','NumberTitle','off','Name','Model');
iconify(1);
subplot(2,1,1);
imagesc(log10(abs(TL.data))');colorbar;
xlabel('datum'),ylabel('time'),title('data(log10)')
refdata=TL.data(:,1);deltadata=TL.data;
for i=1:size(TL.data,2), deltadata(:,i)=TL.data(:,i)./refdata-1; end
subplot(2,1,2);
imagesc(deltadata'*100);colorbar;
xlabel('datum'),ylabel('time'),title('relative deviation in %')
% end

% --------------------------------------------------------------------
function timeinv_Callback(hObject, eventdata, handles)
global Pro Mesh N progdir S TL
refdata=TL.data(:,1);data=TL.data;
for i=1:size(TL.data,2), data(:,i)=log(abs(TL.data(:,i)))-log(refdata); end

err=N.err;
err(min(TL.data)<0)=1000; % downweight wrong signs
D=spdiags(1./log(err+1),0,length(err),length(err));
L=meshsmoothness(Mesh);
DS=D*S;
lambda=str2num(get(handles.lamtext,'String'));
TL.deltamodel=(DS'*DS+30*lambda*L)\(DS'*(D*data));
dc2dtree('timepic_Callback',gcbo,[],guidata(gcbo));        
% end

% --------------------------------------------------------------------
function timepic_Callback(hObject, eventdata, handles)
global Pro Mesh N progdir S TL
figure(1);
set(1,'Units','Characters','MenuBar','none','NumberTitle','off','Name','Model');
iconify(1);set(1,'Color',[1 1 1]);
if isfield(TL,'caxis')&&(length(TL.caxis)>1)&&(diff(TL.caxis)>0),
    cax=TL.caxis;
else
    cax=interperc((exp(TL.deltamodel)-1)*100,[5 95]);
end
if prod(sign(cax))<0,
    cax=[-1 1]*max(abs(cax));cmap=b2r;
else cmap=jet; end
for i=1:size(TL.deltamodel,2),
   clf;colormap(cmap);
   patch('Vertices',Mesh.node,'Faces',Mesh.cell,'FaceVertexCData',...
       (exp(TL.deltamodel(:,i))-1)*100,'FaceColor','flat');%faces,'EdgeColor', edges );
   caxis(cax);
   colorbar;axis equal tight
   title(strrep(TL.names{i},'_',' '));   
   pause(0.1);
   ns=num2str(i);while length(ns)<3, ns=['0' ns]; end
   exportpng1(1,fullfile(TL.dir,[ns '.png']));
%    M(i)=getframe(gcf);M(i).colormap=colormap;
end
% movie2avi(M,fullfile(TL.dir,'time.avi'),'fps',2,'Compression','Indeo3','Quality',100);
% end

% --------------------------------------------------------------------
function timeexport_Callback(hObject, eventdata, handles)
global TL Mesh Pro
if isfield(TL,'first'), 
    [pp,ff,ee]=fileparts(TL.first);
    outfile=fullfile(TL.dir,strrep(TL.first,ee,'.txt'));
else outfile='timelapse.txt'; end
[ff,pp]=uiputfile('*.txt','Select File',outfile);
if ~ischar(ff), return; end
outfile=fullfile(pp,ff);
A=[Mesh.node(Mesh.cell(:,1),1:2),Mesh.node(Mesh.cell(:,3),1:2),Mesh.node(Mesh.cell(:,3),1:2) exp(TL.deltamodel)*100-100];
ss='%.3f';
for i=2:size(A,2), ss=[ss '\t%.3f']; end
ss=[ss '\r\n'];
fid=fopen(outfile,'w');
fprintf(fid,ss,A');
fclose(fid);
% end

% --------------------------------------------------------------------
function timeall_Callback(hObject, eventdata, handles)
global N TL Pro Mesh S
sf=strfind(Pro.dirname,filesep);
alltypes='*.ddz;*.ohm;*.dat';
[TL.first,TL.dir]=uigetfile(alltypes,'select reference file',fullfile(Pro.dirname(1:sf(end)-1),alltypes));
if ~ischar(TL.first), return; end
[pp,ff,ee]=fileparts(TL.first);isddz=0;
fprintf('vor laden');
switch ee,
case '.ddz',
    N0=loadddz(fullfile(TL.dir,TL.first));isddz=1;
otherwise
    N0=sort2delecs(readinv2dfile(fullfile(TL.dir,TL.first),1));
end
if ~isfield(N0,'r'),
    if ~isfield(N0,'rho'),
        if isfield(N0,'u')&isfield(N0,'i'), rho0=N0.u./N0.i; end    
    end
    if isfield(N0,'rho')&isfield(N,'k'), N0.r=N0.rho.*N.k; end    
end
aa=dir(fullfile(TL.dir,['*' ee]));
if isempty(aa), uiwait(errordlg('Did not find files!'));return; end
j=0;TL.names={};TL.data=zeros(length(N.a),1);
abmn0=[N0.a N0.b N0.m N0.n];
L=meshsmoothness(Mesh);
D=spdiags(1./log(N.err+1),0,length(N.err),length(N.err));
% [size(S) size(D)]
DS=D*S;
TL.deltamodel=zeros(Mesh.ncells,1);
for i=1:length(aa),
%     fprintf('Measurement %d\n',i);
    if isddz, NN=loadddz(fullfile(TL.dir,aa(i).name));
    else NN=sort2delecs(readinv2dfile(fullfile(TL.dir,aa(i).name),1)); end
%     NN
    j=j+1;
    if length(aa(j).name)>15, TL.names{i}=aa(j).name(5:15);
    else TL.names{i}=strrep(aa(j).name,ee,''); end
    abmn=[NN.a NN.b NN.m NN.n];
    [c,i0,ii]=intersect(abmn0,abmn,'rows');
    if isfield(NN,'rho'), rho=NN.rho; else rho=NN.u./NN.i; end
%     [size(rho) size(N.k) size(N0.r) length(ii) length(i0)]
    deltadata=log(rho(ii).*N.k(i0))-log(N0.r(i0));
    TL.deltamodel(:,j)=(DS(i0,:)'*DS(i0,:)+30*Pro.lambda*L)\(DS(i0,:)'*(D(i0,i0)*deltadata));
end
cc=inputdlg({'minimum color in % difference (0=automatic)','maximum color in % difference (0=automatic)'},'Specify color axis',1,{'0','0'});
TL.caxis=[str2num(cc{1}) str2num(cc{2})];
dc2dtree('timepic_Callback',gcbo,[],guidata(gcbo));
% end

% --------------------------------------------------------------------
function prepdefault_Callback(hObject, eventdata, handles)
dos('mkdir tmp');
dos('dc2dtreepre -v -B -R 0.15 24.dat');
% end


% --------------------------------------------------------------------
function probeline_Callback(hObject, eventdata, handles)
global Pro Mesh
h=gline;
waitforbuttonpress;waitforbuttonpress;
x=get(h,'XData');
y=get(h,'YData');
phi=atan2(diff(y),diff(x));
ma=[cos(phi) -sin(phi);sin(phi) cos(phi)];
MM=Mesh;MM.node=Mesh.node*ma;
xy=ma'*[x;y];
% line(xy(1,:),xy(2,:))
cellmids=zeros(MM.ncells,2);
for i=1:MM.ncells, cellmids(i,:)=mean(MM.node(MM.cell(i,:),:)); end
xline=linspace(xy(1,1),xy(1,2),30);
rho=griddata1(cellmids(:,1),cellmids(:,2),Pro.res,xline,xy(2,1));
f=figure(1);
set(f,'Units','Pixels','MenuBar','none','NumberTitle','off','Name','Resistivity distribution along line');
iconify(f);
semilogy(xline,rho);grid on;
xlabel('distance in m');ylabel('resistivity in Ohm.m');
set(gca,'YTickLabel',num2strcell(round(get(gca,'YTick'))));

% --------------------------------------------------------------------
function timepickcell_Callback(hObject, eventdata, handles)
global Mesh Pro TL N malstat
if malstat~=0,
    [cmin,cmax]=tripatchmod(Mesh,Pro.res,N);
    malstat=0;
end
cellmids=zeros(Mesh.ncells,2);
for i=1:Mesh.ncells, cellmids(i,:)=mean(Mesh.node(Mesh.cell(i,:),:)); end
[x,y]=ginput(1);
di=sum((cellmids-repmat([x y],Mesh.ncells,1)).^2,2);
[midi,nn]=min(di);
aa=sprintf('Cell No. %d: x=%.3fm y=%.3fm',nn,cellmids(nn,1),cellmids(nn,2));
f=figure(1);clf;
set(f,'Units','Pixels','MenuBar','none','NumberTitle','off','Name','Resistivity variation');
iconify(f);
if isfield(Pro,'res')&&(length(Pro.res)==Mesh.ncells),
    aa=sprintf('%s rho=%.1f Ohmm',aa,Pro.res(nn)); 
    model=(TL.deltamodel(nn,:)+1)*Pro.res(nn);
    plot(TL.dates,model);grid on
    xlabel('time in hours');
    ylabel('\rho in \Omega  m');
    title(aa);
    xt=get(gca,'XTick');
    for i=1:length(xt),
        ss=datestr(xt(i));
        xtl{i}=ss(end-7:end-3);
    end
    set(gca,'XTIckLabel',xtl);
end
% if strcmp(questdlg('Save ),'Yes'),
%     
% end

% if isfield(Pro,'phase')&&(length(Pro.phase)==Mesh.ncells),
%     aa=sprintf('%s phi=%.1f mrad',aa,Pro.phase(nn)); 
% end
% nnn=Mesh.cell(nn,:);nnn(end+1)=nnn(1);
% hold on;line(Mesh.node(nnn,1),Mesh.node(nnn,2),'Color','white');hold off
% iconify(msgbox(aa,'Picked cell'));


% --------------------------------------------------------------------
function batchinversion_Callback(hObject, eventdata, handles)
global N Pro progdir malstat ana S Mesh Invopts Meshopts
set(handles.message,'String',{'DC2dTree - Tree Tomography';'-----------------------------------'});
alltypes='*.ddz;*.ohm;*.dat';
if isfield(Pro,'dirname')&&~isempty(Pro.dirname), 
    ss=strfind(Pro.dirname,filesep);
    infile=fullfile(Pro.dirname(1:ss(end)-1),alltypes);
else
    infile=alltypes;
end
if isfield(Pro,'dirname'), infile=[Pro.dirname filesep '..' filesep alltypes]; end
[fname,pname]=uigetfile({alltypes,'Known files';
    '*.ddz','DDZ Files(*.ddz)';'*.ohm','OHM Files(*.ohm)';...
        '*.dat','DAT Files(*.dat)';'*.*','All files(*.*)'},...
    'Load data file',infile);
if isequal(fname,0), return; end
[pp,ff,ee]=fileparts(fname);
if ~exist(fullfile('tmp','sensMat.mat'),'file')||~exist(fullfile('tmp','meshPara.bms'),'file'),
    dc2dtree('prepdefault_Callback',gcbo,[],guidata(gcbo))
end
Mesh=loadmesh('tmp\meshPara.bms');
MeshP=loadmesh('tmp\meshPrim.bms');
MeshS=loadmesh('tmp\meshSec.bms');
Spp=loadsens;Npp=readinv2dfile('24.data',1);
for i=1:size(Spp,1), Spp(i,:)=Spp(i,:)/Npp.k(i); end
am=[Npp.a(:) Npp.m(:)];    
[midP,angP,oriP]=getdatpar(Npp);
lambda=str2double(get(handles.lamtext,'String'));
%%
set(gcf,'Pointer','watch');
f=figure(1);clf;iconify(f);
set(f,'Units','Pixels','MenuBar','none','NumberTitle','off','Name','Resistivity model');
po=get(1,'Position');po(3)=po(4);set(1,'Position',po);
if strcmp(ee,'.ddz'),
    ss=questdlg('Specify direction of electrode arrangement','Direction','Clockwise','Counterclockwise','Counterclockwise');
    if isequal(ss,'Counterclockwise'), di=-1; else di=1; end    
end
cd(pname);
dos('mkdir tmp');
all=dir(['*' ee]);
for i=1:length(all),
    filename=all(i).name;
    display(['Reading file ' filename '.']);
    switch ee,
        case '.ddz',
            [N,radius]=loadddz(filename,di);
            N.k=kfactorana(N,ana);
            radfile=strrep(filename,'.ddz','rad');
            if exist(radfile,'file'),
                fid=fopen(radfile,'r');
                radius=fscanf(fid,'%f');
                fclose(fid);
                if max(radius)>5, radius=radius/100; end             
                nel=size(N.elec,1);
                while length(radius)<nel, radius(end+1)=radius(end); end
                for i=0:nel-1, %Süden=1, dann Uhrzeigersinn
                    N.elec(i+1,1)=-radius(i+1)*sin(i/nel*2*pi)*di;
                    N.elec(i+1,2)=-radius(i+1)*cos(i/nel*2*pi);
                end            
            end
        case '.dat',
            N=sort2delecs(readinv2dfile(filename,1));
            N.k=kfactorana(N,ana);
            if ~isfield(N,'rho'), N.rho=N.r; end
            N.r=abs(N.rho.*N.k);
        case '.ohm',
            N=sort2delecs(readinv2dfile(filename,1));
            N.k=kfactorana(N,ana);
            if ~isfield(N,'rho'), N.rho=N.r; end
            N.r=abs(N.rho.*N.k);
        end % switch
    rads=sqrt(sum((N.elec-repmat(mean(N.elec),size(N.elec,1),1)).^2,2));
    if (size(N.elec,1)==24)&&(length(unique(round(rads*500)))==1),        
        S=zeros(length(N.a),size(Spp,2));
        for i=1:length(N.a),
            S(i,:)=Spp(find(ismember(am,sort([N.a(i) N.m(i)]),'rows')),:)-...
                Spp(find(ismember(am,sort([N.a(i) N.n(i)]),'rows')),:)-...
                Spp(find(ismember(am,sort([N.b(i) N.m(i)]),'rows')),:)+...
                Spp(find(ismember(am,sort([N.b(i) N.n(i)]),'rows')),:);
            S(i,:)=S(i,:)*N.k(i);
        end
        radius=mean(rads);
        [mid,ang,ori]=getdatpar(N); %mid should be [0 0]
        phi=(ang-angP)/size(N.elec,1)*2*pi;
        twister=[cos(phi) -ori*oriP*sin(phi);ori*oriP*sin(phi) cos(phi)];
        Mesh1=Mesh;Mesh1.node=Mesh.node*radius*twister;
        MeshS1=MeshS;MeshS1.node=MeshS.node*radius*twister;
        MeshP1=MeshP;MeshP1.node=MeshP.node*radius*twister;
        cd(pname);
        savemesh(Mesh1,'tmp\meshPara.bms');
        savemesh(MeshP1,'tmp\meshPrim.bms');
        savemesh(MeshS1,'tmp\meshSec.bms');
        savesensmat(S,'tmp\sensMat.mat');
        cd(progdir);
        datafile='tmp.data';
        D=N;if isfield(D,'rho'), D=rmfield(D,'rho'); end
        if isfield(D,'topo'), D=rmfield(D,'topo'); end
        D.konf=D.k;saveinv2dfile(fullfile(pname,datafile),D);     
    else % standard circle prepare
        cmdline=[progdir filesep 'dc2dtreepre -v '];
        if ~isfield(Meshopts,'nospline')||(Meshopts.nospline==0), cmdline=[cmdline '-B ']; end
        if isfield(Meshopts,'equirefine')&&(Meshopts.equirefine>0), cmdline=[cmdline '-E ']; end
        if isfield(Meshopts,'parasmooth')&&(Meshopts.parasmooth>0), cmdline=[cmdline '-S ']; end
        if isfield(Meshopts,'paraquality'), cmdline=[cmdline '-Q' num2str(Meshopts.paraquality) ' ']; end
        if isfield(Meshopts,'pararefine'), cmdline=[cmdline '-R' num2str(Meshopts.pararefine) ' ']; end
        if isfield(Meshopts,'paramaxarea'), cmdline=[cmdline '-A' num2str(Meshopts.paramaxarea) ' ']; end
        if isfield(Meshopts,'primsmooth')&&(Meshopts.primsmooth>0), cmdline=[cmdline '-s ']; end
        if isfield(Meshopts,'primquality'), cmdline=[cmdline '-q' num2str(Meshopts.primquality) ' ']; end
        if isfield(Meshopts,'primrefine'), cmdline=[cmdline '-r' num2str(Meshopts.primrefine) ' ']; end
        if isfield(Meshopts,'primmaxarea'), cmdline=[cmdline '-a' num2str(Meshopts.primmaxarea) ' ']; end
        if isfield(Meshopts,'secmeshrefine'), cmdline=[cmdline '-f' num2str(Meshopts.secmeshrefine) ' ']; end
        cmdline=[cmdline '"' Pro.datfile '"'];
        set(gcf,'Pointer','watch');
        systemcall(cmdline);        
        Mesh=loadmesh('tmp\meshPara.bms');
    end % individual prepare
    cmdline=[progdir filesep 'dc2dtreerun -v '];
    if get(handles.optlambda,'Value')>0,
        cmdline=[cmdline '-O '];
    else
        cmdline=[cmdline '-l' num2str(lambda) ' '];
    end
    if isfield(Invopts,'lowerbound'), cmdline=[cmdline '-b' num2str(Invopts.lowerbound) ' ']; end
    if isfield(Invopts,'upperbound'), cmdline=[cmdline '-y' num2str(Invopts.upperbound) ' ']; end
    if get(handles.robust,'Value')>0, cmdline=[cmdline '-R ']; end
    if get(handles.blockymodel,'Value')>0, cmdline=[cmdline '-B ']; end
    cmdline=[cmdline '"' datafile '"'];
    cd(pname);
    systemcall(cmdline);
    if exist('tmp\chi2.vector','file'),
        fid=fopen('tmp\chi2.vector');chiq=fscanf(fid,'%f');fclose(fid);
    end
    aa=dir('tmp\model_iter.*');
    cd(progdir);
    for i=1:length(aa), ddd(i)=datenum(aa(i).date); end
    [so,idx]=sort(ddd);
    display(['Taking ' aa(idx(end)).name]);
    fid=fopen(fullfile(pname,'tmp',aa(idx(end)).name));
    Pro.res=fscanf(fid,'%f');fclose(fid);
    figure(1);clf;
    tripatchmod(Mesh1,Pro.res,N);
    epsprint(1,strrep(fullfile(pname,filename),ee,''),1);
    cd(pname);
end % all files
cd(progdir);
Mesh=Mesh1;
dc2dtree('showmod_Callback',gcbo,[],guidata(gcbo))
set(gcf,'Pointer','arrow');
