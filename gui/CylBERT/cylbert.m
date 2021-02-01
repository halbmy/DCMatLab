function varargout = cylbert(varargin)
% CYLBERT M-file for cylbert.fig
%      CYLBERT, by itself, creates a new CYLBERT or raises the existing
%      singleton*.
%
%      H = CYLBERT returns the handle to a new CYLBERT or the handle to
%      the existing singleton*.
%
%      CYLBERT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CYLBERT.M with the given input arguments.
%
%      CYLBERT('Property','Value',...) creates a new CYLBERT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cylbert_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cylbert_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cylbert

% Last Modified by GUIDE v2.5 08-Feb-2008 14:26:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cylbert_OpeningFcn, ...
                   'gui_OutputFcn',  @cylbert_OutputFcn, ...
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


% --- Executes just before cylbert is made visible.
function cylbert_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for cylbert
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(handles.text13,'String','CylBERT 0.7');
set(handles.text23,'String','zweight');
set(handles.zpower,'String',0.01);
global DD dirname
set(handles.columntype,'Value',4);
% dirname='980i';
% dirname='BS335';
ctypes=get(handles.columntype,'String');
ctypes{end+1}='982';
coltype=7;
set(handles.columntype,'String',ctypes);
set(handles.columntype,'Value',coltype);
dirname=ctypes{coltype};
DD=readunifile(fullfile(dirname,'data.def'));

% UIWAIT makes cylbert wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cylbert_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function columntype_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in columntype.
function columntype_Callback(hObject, eventdata, handles)
global DD dirname Data
ss=get(hObject,'String');
dirname=ss{get(hObject,'Value')}
DD=readunifile(fullfile(dirname,'data.def'));
Data=[];
set(handles.file0,'String','');
set(handles.file1,'String','');
set(handles.file2,'String','');
set(handles.file3,'String','');

% --- Executes during object creation, after setting all properties.
function file0_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on button press in dir0.
function dir0_Callback(hObject, eventdata, handles)
global datfile Data
infile='*.txt';
if ~isempty(datfile), 
    [pp,nn,ee]=fileparts(datfile);
    infile=fullfile(pp,'*.txt');    
end
[fname,pname]=uigetfile('*.txt','Read reference data file',infile);
if ~isstr(fname), return; end
pname=strrep(pname,pwd,'.');
datfile=fullfile(pname,fname);
Data0=loadradiccyl(datfile);
if iscell(Data0), %apparently several cycles
    Data=Data0;
    Data0=Data0{1};
%     set(handles.file1,'String',datfile);
%     set(handles.status1,'String',sprintf('%d cycles',length(Data)-1));
else
    Data={Data0};
end    
if ~isstruct(Data0)||(~isfield(Data0,'elec'))||(~isfield(Data0,'a')),
    set(handles.status0,'String','Invalid file');
    return;
end
set(handles.file0,'String',datfile);
ss=sprintf('%d data using %d electrodes',length(Data0.a),size(Data0.elec,1));
if length(Data)>1, 
    ss=[ss ', ' num2str(length(Data)) ' cycles']; 
end
set(handles.status0,'String',ss);

function file0_Callback(hObject, eventdata, handles)

% --- Executes on button press in dir1.
function dir1_Callback(hObject, eventdata, handles)
global datfile Data
[pp,ff,ee]=fileparts(datfile);
[fname,pname]=uigetfile('*.txt','Read cycle data file',fullfile(pp,'*.txt'));
if ~ischar(fname), return; end
pname=strrep(pname,pwd,'.');
datfile=fullfile(pname,fname);
Data1=loadradiccyl(datfile);
%% hotfix
if ~iscell(Data1)&&isstruct(Data1), Data1={Data1}; end
if ~iscell(Data1)||~isstruct(Data1{1})||~isfield(Data1{1},'elec')||~isfield(Data{1},'a'),
    set(handles.status1,'String','Invalid file');
    return;
end
set(handles.file1,'String',datfile);
set(handles.status1,'String',sprintf('%d cycles',length(Data1)));
l=1; %length(Data)
for i=1:length(Data1), 
    if isfield(Data1{i},'elec')&&(size(Data1{i}.elec,1)==size(Data{i}.elec,1))&&...
            isfield(Data1{i},'a')&&(length(Data1{i}.a)==length(Data{1}.a)),
        l=l+1;
        Data{l}=Data1{i};
    end
end
        
% --- Executes during object creation, after setting all properties.
function file1_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function file1_Callback(hObject, eventdata, handles)


% --- Executes on button press in dir2.
function dir2_Callback(hObject, eventdata, handles)
global datfile Data
[pp,ff,ee]=fileparts(datfile);
[fname,pname]=uigetfile('*.txt','Read cycle data file',fullfile(pp,'*.txt'));
if ~ischar(fname), return; end
pname=strrep(pname,pwd,'.');
datfile=fullfile(pname,fname);
Data1=loadradiccyl(datfile);
if ~iscell(Data1)&&isstruct(Data1), Data1={Data1}; end
if ~iscell(Data1)||~isstruct(Data1{1})||~isfield(Data1{1},'elec')||~isfield(Data{1},'a'),
    set(handles.status1,'String','Invalid file');
    return;
end
set(handles.file2,'String',datfile);
set(handles.status2,'String',sprintf('%d cycles',length(Data1)));
l=length(Data);
for i=1:length(Data1), 
    if isfield(Data1{i},'elec')&&(size(Data1{i}.elec,1)==size(Data{i}.elec,1))&&...
            isfield(Data1{i},'a')&&(length(Data1{i}.a)==length(Data{1}.a)),
        l=l+1;
        Data{l}=Data1{i};
    end
end


% --- Executes during object creation, after setting all properties.
function file2_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function file2_Callback(hObject, eventdata, handles)


% --- Executes on button press in dir3.
function dir3_Callback(hObject, eventdata, handles)
global datfile Data
[pp,ff,ee]=fileparts(datfile);
[fname,pname]=uigetfile('*.txt','Read cycle data file',fullfile(pp,'*.txt'));
if ~ischar(fname), return; end
pname=strrep(pname,pwd,'.');
datfile=fullfile(pname,fname);
Data1=loadradiccyl(datfile);
if ~iscell(Data1)&&isstruct(Data1), Data1={Data1}; end
if ~iscell(Data1)||~isstruct(Data1{1})||~isfield(Data1{1},'elec')||~isfield(Data{1},'a'),
    set(handles.status1,'String','Invalid file');
    return;
end
set(handles.file3,'String',datfile);
set(handles.status3,'String',sprintf('%d cycles',length(Data1)));
l=length(Data);
for i=1:length(Data1), 
    if isfield(Data1{i},'elec')&&(size(Data1{i}.elec,1)==size(Data{i}.elec,1))&&...
            isfield(Data1{i},'a')&&(length(Data1{i}.a)==length(Data{1}.a)),
        l=l+1;
        Data{l}=Data1{i};
    end
end

% --- Executes during object creation, after setting all properties.
function file3_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function file3_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function rhoamin_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function rhoamin_Callback(hObject, eventdata, handles)
if ~isfinite(str2double(get(hObject,'String'))), set(hObject,'String','0'); end


% --- Executes during object creation, after setting all properties.
function rhoamax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function rhoamax_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String'))), set(hObject,'String','Inf'); end


% --- Executes during object creation, after setting all properties.
function errormax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function errormax_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String'))), set(hObject,'String','1'); end


% --- Executes during object creation, after setting all properties.
function ipmax_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function ipmax_Callback(hObject, eventdata, handles)
if isnan(str2double(get(hObject,'String'))), set(hObject,'String','10'); end

% --- Executes during object creation, after setting all properties.
function errperc_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function errperc_Callback(hObject, eventdata, handles)
if ~isfinite(str2double(get(hObject,'String'))), set(hObject,'String','2'); end


% --- Executes during object creation, after setting all properties.
function errvolt_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function errvolt_Callback(hObject, eventdata, handles)
if ~isfinite(str2double(get(hObject,'String'))), set(hObject,'String','10'); end


% --- Executes during object creation, after setting all properties.
function defcurr_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function defcurr_Callback(hObject, eventdata, handles)
if ~isfinite(str2double(get(hObject,'String'))), set(hObject,'String','50'); end


% --- Executes on button press in checkerr.
function checkerr_Callback(hObject, eventdata, handles)
global DD Data
errperc=str2num(get(handles.errperc,'String'));
errvolt=str2num(get(handles.errvolt,'String'))*1e-6;
defcurrent=str2num(get(handles.defcurr,'String'))*1e-3;
N=Data{1};
% N.r=N.rho.*DD.k;
err=estimateerror(N,errperc,errvolt,defcurrent);
set(handles.errtext,'String',sprintf('%.1f-%.1f%%',min(err)*100,max(err)*100));

% --- Executes on button press in writefiles.
function writefiles_Callback(hObject, eventdata, handles)
global DD dirname Data
% minrhoa=0;maxerr=0.01;maxip=10;maxrhoa=1
% errperc=2;errvolt=10e-6;defcurrent=50e-3;
errperc=str2num(get(handles.errperc,'String'));
errvolt=str2num(get(handles.errvolt,'String'))*1e-6;
defcurrent=str2num(get(handles.defcurr,'String'))*1e-3;
minrhoa=str2num(get(handles.rhoamin,'String'));
maxrhoa=str2num(get(handles.rhoamax,'String'));
maxerr=str2num(get(handles.errormax,'String'));
maxip=str2num(get(handles.ipmax,'String'));
aa=dir(fullfile(dirname,'*.dat'));
for i=1:length(aa), delete(fullfile(dirname,aa(i).name)); end
if length(Data)>1, fid=fopen(fullfile(dirname,'timesteps.txt'),'w'); end
invalid=[];
abmn=[DD.a(:) DD.b(:) DD.m(:) DD.n(:)];
for i=1:length(Data),    
    fnames{i}=[num2str(i-1) '.dat'];
    N=Data{i};
    N.elec=DD.elec;
%     N.k=DD.k;
%     N.r=N.rho.*N.konf; % alt
    %% ausw�hlen welche (wie timelapse)!!!
    [tf,loc]=ismember([N.a(:) N.b(:) N.m(:) N.n(:)],abmn,'rows');
    fo=find(loc);fo0=find(tf==0);
    N.r=zeros(size(N.a));
    N.k=DD.k(loc(fo));
    N.r(fo)=N.rho(fo).*N.k;    
    N.ip=abs(N.ip)*pi/180*1000;
%     minmax(N.r)
    fi=find((N.r<minrhoa)|(N.r>maxrhoa)|(N.err>maxerr/100)|(N.ip>maxip));
    invalid(i)=length(fi)+length(fo0);
    N.r(fi)=median(N.r);
%     N.rho(fi)=N.r(fi)./N.konf(fi);
    N.err=estimateerror(N,errperc,errvolt,defcurrent);
    N.err(fi)=1e6;N.ip(fi)=0;    
    N.err(fo0)=1e6;N.r(fo0)=median(N.r);N.ip(fo0)=0;
    saveunifile(fullfile(dirname,fnames{i}),N);
%     saveinv3dfile(fullfile(dirname,fnames{i}),rmfield(N,{'r'}));
    if i>1, fprintf(fid,'%s\n',fnames{i}); end
end
if length(Data)>1, fclose(fid); end
set(handles.invalid,'String',sprintf('%d ',invalid));

% --- Executes during object creation, after setting all properties.
function lambda_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function lambda_Callback(hObject, eventdata, handles)
if ~isfinite(str2double(get(hObject,'String'))), set(hObject,'String','20'); end


% --- Executes during object creation, after setting all properties.
function zpower_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function zpower_Callback(hObject, eventdata, handles)
if ~isfinite(str2double(get(hObject,'String'))), set(hObject,'String','0.0'); end

% --- Executes on button press in robust.
function robust_Callback(hObject, eventdata, handles)

% --- Executes on button press in blocky.
function blocky_Callback(hObject, eventdata, handles)

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
global dirname datfile Data
oldpwd=pwd;
try,
    [pp,ff,ee]=fileparts(datfile);
    cd(dirname);
    set(gcf,'Pointer','watch'); 
    if get(handles.lambdaopt,'Value')>0,
        opts='-O ';
    else
        opts=sprintf('-l %s ',get(handles.lambda,'String'));
    end    
    opts=sprintf('%s-z %s ',opts,get(handles.zpower,'String'));
    opts=sprintf('%s-b %s ',opts,get(handles.rhoamin,'String'));
    opts=sprintf('%s-u %s ',opts,get(handles.rhoamax,'String'));
    if get(handles.robust,'Value'), opts=[opts '-R ']; end
    if get(handles.blocky,'Value'), opts=[opts '-B ']; end
    if length(Data)>1, opts=[opts '-t timesteps.txt ']; end
    % wieder einklammern und tauschen
    if exist(fullfile(dirname,'dcinv.result.vtk'),'file'), delete('dcinv.result.vtk'); end
    if exist('primpot.bmat','file'),
        cmdline=['..\dcinv -vST -y primpot.bmat -p meshPara.bms ' opts '0.dat']
    else
        cmdline=['..\dcinv -v -p meshPara.bms ' opts '0.dat']
    end
    dos(cmdline);
    set(gcf,'Pointer','arrow'); 
    if ~exist('dcinv.result.vtk','file'), fprintf('no result!');return; end
    cd(oldpwd);
    if ~exist(fullfile(dirname,'dcinv.result.vtk'),'file'), fprintf('no result!');return; end
    %     dos(['copy /y ' fullfile(dirname,'dcinv.result.vtk') ' "' fullfile(pp,'0.vtk') '"']);
    copyfile(fullfile(dirname,'dcinv.result.vtk'),fullfile(pp,'result.vtk'));
    if length(Data)>1, % timesteps data
        ABS=loadsensmat(fullfile(dirname,'modelAbs.bmat'));
        % DIF=loadsensmat(fullfile(dirname,'modelDiff.mat'));
        Mesh=loadvtkmesh(fullfile(dirname,'dcinv.result.vtk'),3);
        out=sprintf('Model 0: min/max=%.1f/%.1f',min(ABS(1,:)),max(ABS(1,:)));    
        set(handles.invalid,'String',out);
        for i=2:size(ABS),
            di=(ABS(i,:)./ABS(1,:)-1)'*100;
            mesh2vtk(fullfile(pp,['diff' num2str(i-1) '.vtk']),Mesh,di,'diff');
        end
        for i=1:size(ABS,1),
            Mesh.cellattr=ABS(i,:)';
            savevtkmesh(Mesh,fullfile(pp,['abs' num2str(i-1) '.vtk']));
        end
    end
catch,
    cd(oldpwd);
end

% --- Executes on button press in display.
function display_Callback(hObject, eventdata, handles)
global dirname
Mesh=loadmesh(fullfile(dirname,'meshParaDomain.bms'),3);
res=load(fullfile(dirname,'resistivity.vec'));
tetdef=[1 2 3;1 2 4;1 3 4;2 3 4];
cmap=colormap;cach=log10(interperc(res));clog=1;
ind=1:length(res);%find(res<200);
clf;
for i=1:length(ind),
    pts=Mesh.node(Mesh.cell(ind(i),1:4),1:3);
    col=[1 0 0];val=res(ind(i));if clog, val=log10(val); end
    cind=1+round((val-cach(1))/(cach(2)-cach(1))*(length(cmap)-1));
    if cind>length(cmap), cind=length(cmap); end
    if cind<1, cind=1; end
    col=cmap(cind,:);
    for j=1:4,
        td=tetdef(j,:);
        patch(pts(td,1),pts(td,2),pts(td,3),col);
    end
end


% --- Executes on button press in exit.
function exit_Callback(hObject, eventdata, handles)
delete(gcbf);




% --- Executes on button press in lambdaopt.
function lambdaopt_Callback(hObject, eventdata, handles)
soo='On';
if get(handles.lambdaopt,'Value')>0, soo='Off'; end
set(handles.lambda,'Enable',soo);
