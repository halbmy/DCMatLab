function varargout = dc2dinvres(varargin)
% DC2DINVRES Application M-file for dc2dinvres.fig
%    FIG = DC2DINVRES launch dc2dinvres GUI.
%    DC2DINVRES('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 09-Sep-2007 20:24:59
global MAL INV FOR output libmmfile
if nargin == 0  % LAUNCH GUI
    
    fig = openfig(mfilename,'reuse');
   
    % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    guidata(fig, handles);
    iconify(fig);
%     global Mod N MAL INV FOR output malstat datfile
    output=handles.message;
    set(handles.message,'String',...
        {'DC2DInvRes v2.12 - Thomas Günther - RESISTIVITY.NET';'-------------------------------------------'});
    set(handles.inv2d,'Name','DC2dInvRes v2.12 - Thomas Günther');
    set(handles.backtoold,'Enable','Off');    
    if exist('setup.mat','file'), load('setup.mat'); end
    if isempty(MAL),
        MAL=struct('cauto',1,'cmap',0,'clog',1,'xdir',0,'cont',[],...
            'elec',0,'nu',0,'nv',0,'style',0,'showgrid',0,'high',1,'alpha',1);
    end
    if isempty(INV),
        INV=struct('redu',0,'mitschicht',0,'method',0,'auto',2,...
            'start',1,'mico',0.4,'lam',30,'lolo',1,'sens',1,'glob',1,...
            'linesearch',1,'blocky',0,'robust',0,'rbzfak',0.3,'lbound',0,...
            'const',1,'weight',1,'spsens',0);
    end
    if isempty(FOR),
        FOR=struct('rand',4,'prolong',4,'zusatz',2,'method',0,'refine',0,...
            'zref',1,'fillup',1);
    end
    set(handles.preview,'Visible','off');
    libmmfile=checklic2;
    if ~isequal(libmmfile,4),
        set(handles.filesave,'Enable','off');
%         set(handles.fileexport,'Enable','off');
        set(handles.exportmodel,'Enable','off');
        uiwait(warndlg({'No valid license file found! Some functions (export) are not available!',...
                'For a valid license file call Generate License Code from the start menu and send it to thomas@resistivity.net'},...
            'No valid license file')); 
    end
    
    if nargout > 0
        varargout{1} = fig;
    end
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    %else
    
    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
            feval(varargin{:}); % FEVAL switchyard
        end
    catch
        disp(lasterr);
    end
end

% --------------------------------------------------------------------
function varargout = message_Callback(h, eventdata, handles, varargin)

% ---------------------- FILE ----------------------------------------
% --------------------------------------------------------------------
function varargout = file_Callback(h, eventdata, handles, varargin)
global N malstat
set(handles.addlb,'Enable',sonoff(isstruct(N)));
set(handles.fileexport,'Enable',sonoff(~isempty(malstat)));


% --------------------------------------------------------------------
function varargout = fileread_Callback(h, eventdata, handles, varargin)
global Mod N FOR MAL RMS S datfile INV CHIQ libmmfile
alltypes='*.dat;*.ohm;*.flw;*.wen;*.dd;*.pp;*.pd;*.slm;*.amp;*.tx0;*.s2d;*.txt;*.stg;*.bin';
if isempty(datfile), infile=alltypes; else
    [pname,fname,ext]=fileparts(datfile);
    infile=strrep(datfile,[fname ext],alltypes);
end
[filename,pathname]=uigetfile({alltypes,'Known Files';...
    '*.dat/ohm','Data files (*.dat/ohm)';'*.amp','ABEM files (*.amp)';...
    '*.txt','RESECS files (*.txt)';'*.flw;*.dd;*.wen;*.pd;*.slm;*.pp','GEOTOM files (*.flw)';...
    '*.tx0','Lippmann files (*.tx0)';'*.bin','Syscal Pro file (*.bin)';...
    '*.stg','Sting device files (*.stg)';'*.s2d','Soundings file (*.s2d)';...
    '*.*','All files'},'Read Data File',infile);
if filename==0, return; end
set(handles.message,'String',...
    {'DC2DInvRes v2.12 - Thomas Günther';'----------------'});
datfile=strrep(fullfile(pathname,filename),[pwd filesep],'');
if ~exist(datfile,'file'),
    error(['Could not read datafile! ' datfile]);
    return; 
end
libmmfile=checklic2;
set(gcf,'Pointer','watch');
N=read2dfile(datfile);
if (size(N.elec,2)>2)&&(max(N.elec(:,3))>eps), 
    % other case (x=const, y=x)?
    N.topo=N.elec(:,[1 3]);N.elec(:,2)=0;N.elec(:,3)=[];
end
if ~isfield(N,'a')||isempty(N.a), 
    errordlg(['Error opening datafile ' datfile],'Wrong File');
    return;
end
fi=find(N.r<0);ndata=length(N.a);
if length(fi)/ndata>0.5, % apparently systematicly wrong sign
    N.r=abs(N.r);N.rho(fi)=-N.rho(fi);
end
fi=find(N.r<0);ndata=length(N.a);
if ~isempty(fi),
    N=delmeasurement(N,fi);
    messg(sprintf('Found negative app. res.! Neglecting %d data',length(fi)));
end
fi=find(isinf(N.r)|isnan(N.r));
if ~isempty(fi),
    N=delmeasurement(N,fi);
    messg(sprintf('Found invalid (Inf,NaN) app. res.! Neglecting %d data',length(fi)));
end
dc2dinvres('showdata_Callback',gcbo,[],guidata(gcbo));
if (~isfield(N,'err'))||isempty(N.err)||(min(N.err)<0),
    daterr;
    dc2dinvres('showdata_Callback',gcbo,[],guidata(gcbo));
else
    messg(sprintf('Found error estimate in file, min=%.1f%%, max=%.1f%%',min(N.err)*100,max(N.err)*100));
    if min(N.err)<0.01,
        aa=inputdlg('Type percentage error to be added','Small errors found (2%)!',1,{'1'});
        bb=str2num(aa{1});
        if isempty(bb), bb=2; end
        N.err=N.err+bb/100;
    end
end
[fi1,fi2]=getrez(N);
if any(fi1),
    messg(sprintf('Found reciprocals for %d data',length(fi1)));
    R1=N.r(fi1);
    R2=N.r(fi2);
    set(figure(9),'MenuBar','none','NumberTitle','off','Name','Reciprocity Crossplot');
    loglog(R1,R2,'.');xlabel('normal');ylabel('reverse');
    rez=(R1-R2)./(R1+R2);
    messg(sprintf('std/max reciprocity = %.1f/%.1f%%',...
        std(rez)*100,max(rez)*100));
    N.rez=N.a*0;
    N.rez(fi1)=rez;
    N.rez(fi2)=rez;
    drawnow;
end
% Create & plot model
[fpath,name,ext]=fileparts(datfile);
matfile=strrep(datfile,ext,'-sens.mat');
if exist(matfile,'file'), % found sens. matrix
    libmmfile=checklic2;
    xsave=[];zsave=[];
    load(matfile);
    if any(xsave)&&any(zsave),         
        Mod.x=round(xsave*1000)/1000;Mod.z=round(zsave*1000)/1000;
        messg('Loading sensitivity and parameterization from file');
    else
        [Mod.x,Mod.z]=modelfromdata2d(N); 
    end
    rq=median(N.r);
    Mod.M=ones(length(Mod.x)-1,length(Mod.z)-1)*rq;
    Mod.Mref=Mod.M;
    Mod.R=ones(size(N.r))*rq;Mod.isfor=1;
    if isfield(N,'ip')&&(length(N.ip)==length(N.a)), Mod.ipfor=zeros(size(N.a)); end
    Mod.Cov=zeros(size(Mod.M));
    nm=numel(Mod.M);
    if size(S,2)==nm,
        %     Mod.Cov(:)=sum(abs(S));
        for i=1:nm, Mod.Cov(i)=sum(abs(S(:,i))); end
    end
    if size(S,1)==nm,
        %     Mod.Cov(:)=sum(abs(S),2);
        for i=1:nm, Mod.Cov(i)=sum(abs(S(i,:))); end
    end
    messg(strcat(sprintf('Model(%dx%d=%d cells): x=%g..%g',size(Mod.M,1),size(Mod.M,2),numel(Mod.M),min(Mod.x),max(Mod.x)),...
        '  z= ',sprintf('%g ',Mod.z),sprintf('  Rho_0 = %.1f',rq)));  
    figure(2); %to enforce comparedata window
    dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
    previewon(handles);
else
    [Mod.x,Mod.z,Mod.M]=modelfromdata2d(N);
    if length(Mod.z)>length(Mod.x), % possibly crosshole case
        dx=median(diff(Mod.z));
        Mod.x=(min(N.elec(:,1))-5*dx:dx:max(N.elec(:,1))+5*dx)';
        rho=Mod.M(1,1);Mod.M=ones(length(Mod.x)-1,length(Mod.z)-1)*rho;
        if isfield(N,'ip')&&(length(N.ip)==length(N.a)), Mod.ipfor=zeros(size(N.a)); end
        Mod.Mref=Mod.M;
    end
    dx=median(diff(sort(N.elec(:,1))));
    while length(Mod.x)>size(N.elec,1)*2, % wrong electrode distance
        dx=dx*2;dx=round(dx);
        xmin=min(N.elec(:,1))-dx;xmin=round(xmin/dx)*dx;
        xmax=max(N.elec(:,1))+dx;xmax=round(xmax/dx)*dx;
        Mod.x=(xmin:dx:xmax)';
    end
%     if length(find(N.elec(:,2)==2))<5, dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo)); end
    messg(strcat(sprintf('Model(%dx%d=%d cells): Mod.x=%g..%g',size(Mod.M,1),size(Mod.M,2),numel(Mod.M),min(Mod.x),max(Mod.x)),...
        '  Mod.z= ',sprintf('%g ',Mod.z),sprintf('  Rho_0 = %.1f',Mod.M(1,1))));  
    S=[];Mod.Cov=[];
    dc2dinvres('loadsens_Callback',gcbo,[],guidata(gcbo)); %???
end
rho=Mod.M(1,1);Mod.Lay=rho;
% Starting model (halfspace)
Mod.M=ones(length(Mod.x)-1,length(Mod.z)-1)*rho;Mod.Mref=Mod.M;
% Forward Calculation number "zero"
Mod.R=ones(size(N.r))*rho;Mod.isfor=1;
RMS=rms(N.r,Mod.R,INV.lolo);
CHIQ=chi2(N.r,Mod.R,N.err,INV.lolo);
dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
dc2dinvres('clearsvd_Callback',gcbo,[],guidata(gcbo))
try 
    set(h,'Name',[ 'DC2DInvRes - ' datfile ]);
catch
    set(handles.inv2d,'Name',[ 'DC2DInvRes - ' datfile ]);
end
if exist('model.mat','file'), delete('model.mat'); end
set(handles.backtoold,'Enable','Off');    
set(gcf,'Pointer','arrow');
%if INV.method==0, INV.lam=0; end

% --------------------------------------------------------------------
function varargout = filewrite_Callback(h, eventdata, handles, varargin)
global datfile N
[fpath,name,ext]=fileparts(datfile);
outfile=[];
if ~isempty(datfile), outfile=strrep(datfile,ext,'.dat'); end
if isempty(outfile), outfile='*.dat'; end
[fname,pname]=uiputfile(outfile,'Save Datum Points');
if fname~=0,
    datfile=strrep(fullfile(pname,fname),[pwd filesep],'');
    [ff,pp,ee]=fileparts(datfile);
    if isempty(ee)||strcmp(ee,''), datfile=[datfile '.dat']; end
    messg(['Saving data to file ',datfile]);
    saveinv2dfile(datfile,N);
    datfile=strrep(datfile,[pwd filesep],'');
    try 
        set(h,'Name',[ 'DC2DInvRes - ' datfile ]);
    catch
        set(handles.inv2d,'Name',[ 'DC2DInvRes - ' datfile ]);
    end
end

% --------------------------------------------------------------------
function varargout = addlb_Callback(h, eventdata, handles, varargin)
global Mod N FOR MAL RMS S datfile INV CHIQ
alltypes='*.dat;*.ohm;*.flw;*.amp;*.tx0;*.s2d;*.txt;*.stg';
if isempty(datfile), infile=alltypes; else
    [pname,fname,ext]=fileparts(datfile);
    infile=strrep(datfile,[fname ext],alltypes);
end
[filename,pathname]=uigetfile({alltypes,'Known Files';...
    '*.dat/ohm','Data files (*.dat/ohm)';'*.amp','ABEM files (*.amp)';...
    '*.txt','RESECS files (*.txt)';'*.flw','GEOTOM files (*.flw)';...
    '*.tx0','Lippmann files (*.tx0)';...
    '*.stg','Sting device files (*.stg)';'*.s2d','Soundings file (*.s2d)';...
    '*.*','All files'},'Read Data File',infile);

% [pname,fname,ext]=fileparts(datfile);
% [filename,pathname]=uigetfile({'*.dat';'*.*'},'Read Data File',strrep(datfile,[fname ext],'*.dat'));
if filename==0, return; end
%[path,oname,ext]=datfile;
datfile=strrep(fullfile(pathname,filename),[pwd filesep],'');
if ~exist(datfile,'file'), error('Could not read datafile!'); end
S=[];Mod.Cov=[];
N1=N;
N=read2dfile(datfile);

if (~isfield(N,'err'))||(min(N.err)<=0),
    daterr;
end
N=combdata2d(N,N1);
messg(sprintf('Sum: %d Electrodes, %d Datapoints',size(N.elec,1),length(N.a)));
[Mod.x,Mod.z,Mod.M]=modelfromdata2d(N);

rho=Mod.M(1,1);Mod.Lay=rho;
% Starting model (halfspace)
Mod.M=ones(length(Mod.x)-1,length(Mod.z)-1)*rho;Mod.Mref=Mod.M;
if isfield(N,'ip')&&(length(N.ip)==length(N.a)), Mod.ipfor=zeros(size(N.a)); end
dc2dinvres('showdata_Callback',gcbo,[],guidata(gcbo))
% Forward Calculation number "zero"
Mod.R=ones(size(N.r))*rho;Mod.isfor=1;
RMS=rms(N.r,Mod.R,INV.lolo);
CHIQ=chi2(N.r,Mod.R,N.err,INV.lolo);
dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
dc2dinvres('clearsvd_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function varargout = savelb_Callback(h, eventdata, handles, varargin)
global datfile N
[path,name,ext]=fileparts(datfile);
outfile=strrep(datfile,ext,['1' ext]);
[fname,pname]=uiputfile(outfile,'Save Datum Points');
if fname~=0,
    datfile=strrep(fullfile(pname,fname),[pwd filesep],'');
    [ff,pp,ee]=fileparts(datfile);
    if strcmp(ee,''), datfile=[datfile '.dat']; end
    messg(['Exporting to Res2dinv file ',datfile]);
    if isfield(N,'ip'), saveres2dinvfile(datfile,N,N.r,N.ip);
    else saveres2dinvfile(datfile,N); end
    try 
        set(h,'Name',[ 'DC2DInvRes - ' datfile ]);
    catch
        set(handles.inv2d,'Name',[ 'DC2DInvRes - ' datfile ]);
    end
end

% --------------------------------------------------------------------
function varargout = filesave_Callback(h, eventdata, handles, varargin)
global Mod N MAL INV FOR RMS datfile S CHIQ FIX XX ZZ
mess=get(handles.message,'String');
[path,name,ext]=fileparts(datfile);
matfile=strrep(datfile,ext,'.mat');
[fname,pname]=uiputfile(matfile,'Save Workspace');
if fname~=0, % Testversion nächste Zeile einklammern
%     errordlg('Save workspace not supported in Test version!');return;
    matfile=fullfile(pname,fname);
    messg(['Saving workspace to file ' matfile]);
    save(matfile);
end

% --------------------------------------------------------------------
function varargout = fileload_Callback(h, eventdata, handles, varargin)
global Mod N MAL INV FOR RMS datfile S CHIQ libmmfile FIX XX ZZ
mess=[];
if isempty(datfile)
    matfile='*.mat';
else
    [pp,nn,ee]=fileparts(datfile);
    matfile=strrep(datfile,ee,'.mat');
end
[fname,pname]=uigetfile(matfile,'Load Workspace');
if fname~=0,
    matfile=fullfile(pname,fname);
    libmmfile=checklic2;
    load(matfile);
    if isempty(Mod), % compatibility to v<2.12
        Mod.x=x;Mod.z=z;Mod.M=M;Mod.Cov=zeros(size(Mod.M));
        Mod.R=R;Mod.Lay=Lay;Mod.Mref=Mod.M; end
    if isempty(datfile),
        datfile=strrep(matfile,'.mat','.dat');
        datfile=strrep(datfile,[pwd filesep],'');
    end
    if (~isempty(mess))&&ishandle(handles.message), 
        set(handles.message,'String',mess); 
    else
        messg('(Loading workspace from disk file)');
        dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
    end
    try
        set(h,'Name',[ 'DC2DInvRes - ' datfile ]);
    catch
        set(handles.inv2d,'Name',[ 'DC2DInvRes - ' datfile ]);
    end
    dc2dinvres('loadsens_Callback',gcbo,[],guidata(gcbo))
    dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
end

% --------------------------------------------------------------------
function varargout = fileexport_Callback(h, eventdata, handles, varargin)
global datfile malstat N Mod MAL S INV VD RM RD rmn libmmfile FIX XX ZZ
f=figure(9);set(9,'Menubar','none');
set(f,'Color',[1 1 1]);
clf;
pcolor(rand(3));
%alpha(rand(2));
cb=colorbar('horiz');
set(gca,'Units','normalized');
set(cb,'Units','normalized');

set(handles.malfeld,'Units','pixel');
poall=get(handles.inv2d,'Position');
rel=get(handles.message,'Position');rel=rel(4)/poall(4); %new!
pof=get(gcbf,'position');pof(4)=pof(4)*(1-rel);
%pof(2)=pof(2)-rel;
set(f,'Units','Character','position',pof); % gleiche figure-groesse
set(handles.malfeld,'Units','normalized'); % zurück
po=get(handles.malfeld,'position');
po1=get(gca,'Position');
ind=[1 3 4];po1(ind)=po(ind);
set(gca,'Position',po1);

[path,name,ext]=fileparts(datfile);
if isempty(malstat), malstat=0; end
if ismember(malstat,[6 7 8 10 11]),
    if isempty(VD),
        dc2dinvres('compsvd_Callback',gcbo,[],guidata(gcbo))
    end
end
if malstat>8,
    mal=MAL;
    mal.log=0;mal.cmap=2;mal.cauto=0;
    nr=round(get(handles.slider,'Value'));
end
switch malstat,
    case 0, % Model
%         Mod.Cov=Mod.M*0;if size(S,2)==prod(size(Mod.M)), Mod.Cov(:)=sum(abs(S)); end
        mal=MAL;mal.canot='Ohm*m';
        if max(Mod.Cov(:))>0, patch2dmodel(Mod.x,Mod.z,Mod.M,mal,N,Mod.Cov);
        else patch2dmodel(Mod.x,Mod.z,Mod.M,mal,N); end
        plotconstraints(Mod,FIX,XX,ZZ,N);    
        zus='';
    case -1, % Error
        mal=MAL;mal.cauto=1;
        showdata2d(N,N.err*100,mal);
        zus='-error';
    case 1, % Data
        showdata2d(N,N.r,MAL);
        zus='-data';
    case 2, % Forward
        showdata2d(N,Mod.R,MAL);
        zus='-forw';
    case 3, % Misfit
        dR=(Mod.R(:)./N.r(:)-1)*100;
        showdata2d(N,dR,MAL);
        zus='-misfit';
    case 4, % Coverage
        ma=struct('clog',1);
%         draw2dmodel(Mod.x,Mod.z,log10(sum(abs(S))),ma);
        patch2dmodel(Mod.x,Mod.z,Mod.Cov,ma,N);
        zus='-cov';
    case 5, % DOI-Index
%         draw2dmodel(Mod.x,Mod.z,getdoi(N.r,S,Mod.M(1,1),1.05,INV));
        patch2dmodel(Mod.x,Mod.z,getdoi(N.r,S,Mod.M(1,1),1.05,INV),[],N);
        zus='-doi';
    case 6, % Model resolution
        Rm=Mod.M;
        Rm(:)=diag(RM);
        mal=MAL;
        mal.cauto=0;mal.log=0;mal.cmin=0;mal.cmax=1;
%         draw2dmodel(Mod.x,Mod.z,Rm,mal);
        patch2dmodel(Mod.x,Mod.z,Rm,mal,N);
        zus='-modres';
    case 7, % Data resolution
        Rd=diag(RD);
        mal=MAL;
        mal.cauto=0;mal.log=0;mal.cmin=0;mal.cmax=1;
        showdata2d(N,Rd,mal);
        zus='-datres';
    case 8, % Resolution radius
        Rm=diag(RM);
        [DX,DZ]=ndgrid(diff(Mod.x),diff(Mod.z));
        Fl=DX.*DZ;
        Resrad=sqrt(Fl(:)./Rm(:)/pi);
        mal=MAL;mal.cauto=1;mal.log=1;
%         draw2dmodel(Mod.x,Mod.z,Resrad,mal);
        patch2dmodel(Mod.x,Mod.z,Resrad,mal,N);
        zus='-resrad';
    case 9, % Sensitivity
        snr=Mod.M;snr(:)=S(nr,:);
%         mm=max(abs(snr))/10;
        mm=std(snr(:))*3;
        mal.cmin=-mm;mal.cmax=mm;
%         draw2dmodel(Mod.x,Mod.z,snr,mal);
        patch2dmodel(Mod.x,Mod.z,snr,mal,N);
        zus=['-sens' num2str(nr)];
    case 10, % Data Covariance
        [mids,seps,ii,kk]=showdata2d(N,RD(:,nr));
        hold on;plot(mids(ii(nr)),kk(nr),'kx','MarkerSize',10); hold off
        zus=['-rd' num2str(nr)];
    case 11, % Model Covariance
        mm=max(abs(rmn));mal.cmin=-mm;mal.cmax=mm;
%         draw2dmodel(Mod.x,Mod.z,rmn,mal);
        patch2dmodel(Mod.x,Mod.z,rmn,mal,N);
        lx=length(Mod.x)-1;
        nx=mod(nr,lx);if nx==0, nx=lx; end
        nz=(nr-nx)/lx+1;
        hold on
        line([Mod.x(nx) Mod.x(nx+1) Mod.x(nx+1) Mod.x(nx) Mod.x(nx)],[Mod.z(nz) Mod.z(nz) Mod.z(nz+1) Mod.z(nz+1) Mod.z(nz)],...
            'LineWidth',2,'Color','Black');
        hold off
        zus=['-rm' num2str(nr)];
    case 12, % Voltage
        mal=MAL;mal.cauto=1;
        showdata2d(N,abs(N.u)*1000,struct('log',1));
        zus='-u';
    case 13, % Current
        mal=MAL;mal.cauto=1;
        showdata2d(N,N.i*1000,mal);
        zus='-i';
    case 14, % ip data
        mal=MAL;mal.cauto=1;mal.clog=0;
        showdata2d(N,N.ip,mal);
        zus='-ipdata';
    case 15, % Chargeability
        mal=MAL;mal.cauto=1;
        mal.clog=(min(Mod.IP(:))>0);mal.log=mal.clog;
        if max(Mod.Cov(:))>0, patch2dmodel(Mod.x,Mod.z,Mod.IP,mal,N,Mod.Cov);
        else patch2dmodel(Mod.x,Mod.z,Mod.IP,mal,N); end
        plotconstraints(Mod,FIX,XX,ZZ,N);    
        zus='-ipmodel';
    otherwise
        close(f);
        return;
end
if (malstat==9)||(malstat==10), % Sens. oder Single Datres
    ch='AM';
    el=[N.elec(N.a(nr),:);N.elec(N.m(nr),:)];
    if N.b(nr)>0, el=[el;N.elec(N.b(nr),:)];ch=[ch 'B']; end
    if N.n(nr)>0, el=[el;N.elec(N.n(nr),:)];ch=[ch 'N']; end
    if malstat==10, el(:,2)=0.5; end
    xt=get(gca,'XTick');xl=cellstr(get(gca,'XTickLabel'));
    for k=1:length(ch),
        l=find(xt==el(k,1));
        if isempty(l), l=length(xt)+1; end
        xl{l}=ch(k);xt(l)=el(k,1);
    end
    [xt,l]=sort(xt);
    yl=xl;for i=1:length(xl), yl{i}=xl{l(i)}; end
    set(handles.malfeld,'XTick',xt,'XTickLabel',yl);
    hold on
    plot(el(:,1),el(:,2),'kv','MarkerSize',2);
    hold off
end

%set(f,'Renderer','painters','RendererMode','manual');
set(f,'Units','pixel','PaperSize',po(3:4),'PaperPositionMode','auto');
set(gca,'YTickMode','manual','YTickLabelMode','manual');
exts={'*.eps ';'*.png ';'*.pdf '};
zus=strcat(zus,exts{1});
outfile=strrep(datfile,ext,zus);
[fname,pname]=uiputfile(exts,'Save Figure as',strrep(outfile,'*',''));
if fname~=0,
    outfile=fullfile(pname,fname);
    [ff,pp,ee]=fileparts(outfile);
    if strcmp(ee,''), outfile=[outfile exts{1}]; end
    messg(['Exporting Figure to image ' outfile]);
    ispdf=any(strfind(fname,'.pdf'));
    iseps=any(strfind(fname,'.eps'))|ispdf;
    if iseps,% Testversion nächste Zeile einklammern
        if ~isequal(libmmfile,4),
            close(f);errordlg('No eps/pdf export in Test version');return;
        end
        if ispdf, outfile=strrep(outfile,'.pdf','.eps'); end
        print(f,'-depsc2','-painters',outfile);
        if ispdf, dos(['epstopdf "' outfile '"']);pause(1.0);delete(outfile); end
    else
        exportpng(f,outfile);
%        write_png(getfield(getframe(f),'cdata'),[],outfile);
%         saveas(f,outfile,'png');%print(f,'-r100','-dpng',outfile);
    end
    if malstat==0, % model->additional messages in text file
        if ishandle(2), % data & response too
            endung='png';if iseps, endung='eps'; end
            coutfile=strrep(outfile,['.' endung],['-comp.' endung]);
            messg(['Exporting Data&Response to image ' coutfile]);
            if iseps, 
                print(2,'-depsc2','-painters',coutfile);
                if ispdf, 
                    dos(['epstopdf "' coutfile '"']);
                    pause(1.0);
                    delete(coutfile); 
                end
            else exportpng(2,coutfile); end
        end
        %fid=fopen(strrep(outfile,'.png','.txt'),'w');
        fid=fopen([outfile '.txt'],'w');
        if fid>=0,
            ss=get(handles.message,'String');
            fprintf(fid,'%s\n',ss{:});
            fprintf(fid,'Inversion Options\n');
            smethod={'Matrix','TSVD','SIRT','TLS'};
            sweight={'Equal','Smooth1','Smooth2-dir','Smooth2-neu','Smooth2-mod','Coverage','resolution'};
            sauto={'manual','lcurve','fixed'};
            sredu={'none','delete','combine'};
            fprintf(fid,'%s inversion,weight=%s,regu=%s,reduction=%s,lambda=%.1f,keep=%s\n',...
                smethod{INV.method+1},sweight{INV.weight+1},sauto{INV.auto+1},...
                sredu{INV.redu+1},INV.lam,sauto{INV.const+1});
            %fprintf(fid,'Forward Options\n');
            fclose(fid);
        end
    end
end
%delete(cb);
close(f);

% --------------------------------------------------------------------
function varargout = fileexit_Callback(h, eventdata, handles, varargin)
global output datfile FOR INV MAL
if exist('setup.mat','file'),
    load('setup.mat','FOR','INV','MAL');%,'datfile');
end
save('setup.mat','FOR','INV','MAL','datfile');
output=[];
if exist('model.mat','file'), delete('model.mat'); end
if exist('rhoeta.mat','file'), delete('rhoeta.mat'); end
for i=1:15,
    if ishandle(i), close(i); end
end
delete(gcbf);

% ----------------------- MODEL --------------------------------------
% --------------------------------------------------------------------
function varargout = start_Callback(h, eventdata, handles, varargin)
% the model main menu
global Mod N
set(handles.widerstand,'String','0');
smod=sonoff(isfield(Mod,'M')&&~isempty(Mod.M));
set(handles.modpar,'Enable',smod);
set(handles.exportmodel,'Enable',smod);
set(handles.comparemod,'Enable',smod);
set(handles.blockify,'Enable',smod);
smod=sonoff(isfield(Mod,'M')&&(~isempty(Mod.M)||isstruct(N)));


% --------------------------------------------------------------------
function varargout = startlay_Callback(h, eventdata, handles, varargin)
global N Mod FOR MAL RMS INV S CHIQ
if ~isstruct(N), errordlg('No data present!');return; end
messg('Estimating 1D-Starting Model...');
rho=Mod.Lay(1);
% Mod.Lay=onedinv(N,Mod.z,Mod.M);
[Mod.Lay,Mod.R]=full1dinv(N,Mod.z);
Mod.M(:)=rho;
for i=1:size(Mod.M,2), Mod.M(:,i)=Mod.Lay(i); end
dM=log(Mod.M(:))-log(rho);
%Mod.Mref=Mod.M; % eigentlich keine gute Idee, Mod.z.B. bei smoothness
% dM=log(Mod.M(:))-log(rho);
messg(sprintf('%.1f ',Mod.Lay));
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
dc2dinvres('forwarding_Callback',gcbo,[],guidata(gcbo))
% set(gcf,'Pointer','watch');drawnow;
% if max(N.elec(:,2))>0, %problem with 1d and subsurface electrodes
%     Mod.M(end,end)=Mod.M(end,end)*1.0001;
%     Mod.R=dcfwd2d(Mod.x,Mod.z,Mod.M,Mod.Lay,N,FOR);
% else
%     if length(unique(Mod.Lay))==1,
%         Mod.R=ones(size(N.r))*Mod.Lay(1);
%     else
%         messg('Calculating forward response of 1D model...');
%         Mod.R=fwd2d1d(N,Mod.Lay,Mod.z);
%     end
% end
% set(gcf,'Pointer','arrow');
% RMS=[RMS rms(N.r,Mod.R,INV.lolo)];
% CHIQ=[CHIQ chi2(N.r,Mod.R,N.err,INV.lolo)];
% dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
% if ishandle(2),
%     dc2dinvres('compare_Callback',gcbo,[],guidata(gcbo))
% end
%     S = S + ((log(Mod.R)-log(rho)-S*dM)*(dM'))/(dM'*dM); % broyden
if (size(S,1)==length(Mod.R))&&(size(S,2)==length(dM)),
    ddr=(log(Mod.R)-log(rho)-S*dM(:))/(dM(:)'*dM(:));    
    for i=1:size(S,2), S(:,i)=S(:,i)+ddr*dM(i); end
    messg('Doing Broyden update of sensitivity matrix!');
end
if (size(S,2)==length(Mod.R))&&(size(S,1)==length(dM)), %trans+sparse
    ddr=(log(Mod.R)-log(rho)-dM(:)'*S)/(dM(:)'*dM(:));    
    for i=1:size(S,1), S(i,:)=S(i,:)+ddr*dM(i); end
    messg('Doing Broyden update of sensitivity matrix!');
end
set(handles.backtoold,'Enable','Off');    
%if INV.method==0, INV.lam=0; end

% --------------------------------------------------------------------
function varargout = starthom_Callback(h, eventdata, handles, varargin)
global Mod N FOR MAL RMS INV datfile S CHIQ malstat FIX
if (~isfield(Mod,'M')||isempty(Mod.M))&&isstruct(N), [Mod.x,Mod.z,Mod.M]=modelfromdata2d(N); end
% minkonf=min(abs(N.k(:)));
% rho=median(N.r(abs(N.k)==minkonf));
if isequal(size(FIX),size(Mod.M)),
    [ii,jj]=find(FIX);
    rho=Mod.M(ii(1),jj(1));
elseif isstruct(N), rho=median(N.r); 
else rho=median(Mod.M(:)); end
if rho>20, rho=round(rho/10)*10; end
Mod.Lay=rho;
% Starting model (halfspace)
Mod.M=ones(length(Mod.x)-1,length(Mod.z)-1)*rho;Mod.Mref=Mod.M;
if isfield(N,'ip')&&(length(N.ip)==length(N.a)), Mod.ipfor=zeros(size(N.a)); end
if isstruct(N),
    % Forward Calculation number "zero"
    Mod.R=ones(size(N.r))*rho;Mod.isfor=1;
    RMS=rms(N.r,Mod.R,INV.lolo);
    CHIQ=chi2(N.r,Mod.R,N.err,INV.lolo);
    dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
    S=[];Mod.Cov=[];
%     [path,name,ext]=fileparts(datfile);
    dc2dinvres('loadsens_Callback',gcbo,[],guidata(gcbo))
    %set(handles.message,'String',{'DC2DINVRES v1.0 - Thomas Günther','----------------'});
    %if INV.method==0, INV.lam=0; end
end
INV.first=0;
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
set(handles.backtoold,'Enable','Off');    
previewon(handles);

% --------------------------------------------------------------------
function varargout = setasdata_Callback(h, eventdata, handles, varargin)
global N Mod
N.r=Mod.R;
if isfield(Mod,'ipfor')&&(length(Mod.ipfor)==length(N.a)), N.ip=Mod.ipfor; end
dc2dinvres('showdata_Callback',gcbo,[],guidata(gcbo))

% ----------------- OPTIONS ------------------------------------------
% --------------------------------------------------------------------
function varargout = opt_Callback(h, eventdata, handles, varargin)

% --------------------------------------------------------------------
function varargout = optinv_Callback(h, eventdata, handles, varargin)
s_inv;

% --------------------------------------------------------------------
function varargout = optfor_Callback(h, eventdata, handles, varargin)
s_for;

% --------------------------------------------------------------------
function varargout = optgraf_Callback(h, eventdata, handles, varargin)
fig=s_mal;
uiwait(fig);
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
%figure(fig);

% --------------------------------------------------------------------
function varargout = saveopts_Callback(h, eventdata, handles, varargin)
global FOR MAL INV 
save('setup.mat');

% --------------------------------------------------------------------
function varargout = resetopts_Callback(h, eventdata, handles, varargin)
global FOR MAL INV
MAL=struct('cauto',1,'cmin',100,'cmax',500,'cmap',0,'log',1,'xdir',0,'cont',[],'elec',0);
INV=struct('redu',0,'mitschicht',0,'method',2,'auto',1,...
    'start',1,'mico',0.4,'lam',0.3,'lolo',1,'sens',1);  
FOR=struct('rand',4,'prolong',5,'zusatz',2,'refine',1);
save('setup.mat','FOR','INV','MAL');

% ---------------- INV ----------------------------------------------- 
% --------------------------------------------------------------------
function varargout = inv_Callback(h, eventdata, handles, varargin)
global N Mod
isdata=isstruct(N)&&isfield(N,'a');
isip=isdata&&isfield(N,'ip')&&(length(N.ip)==length(N.a));
set(handles.ipinv,'Enable',sonoff(isip));
ismod=isfield(Mod,'M')&&~isempty(Mod.M);
set(handles.invauto,'Enable',sonoff(isdata&ismod));
set(handles.invone,'Enable',sonoff(isdata&ismod));
set(handles.rollalong,'Enable',sonoff(isdata&ismod));
set(handles.timelapse,'Enable',sonoff(isdata&ismod));
set(handles.startlay,'Enable',sonoff(isdata));


% --------------------------------------------------------------------
function varargout = invauto_Callback(h, eventdata, handles, varargin)
global Mod N FOR MAL RMS S INV CHIQ
%% check lower/upper bound (for broyden update)
lbound=0;ubound=0;
if isfield(INV,'lbound'), lbound=INV.lbound; end
if min(Mod.M(:))<lbound, lbound=0; end
if isfield(INV,'ubound'), ubound=INV.ubound; end
if ubound<max(max(N.r),max(Mod.R)), ubound=0; end
if max(Mod.M(:))>ubound, ubound=0; end
%% run invone until Chi^2<1, deltaChi^2<5% 
running=1;
while running,
    oldMod=Mod;
    dc2dinvres('invone_Callback',gcbo,[],guidata(gcbo))    
    dchiq=CHIQ(end)-CHIQ(end-1);
    if (dchiq<0)||(length(RMS)<3),
        running=(abs(dchiq)/CHIQ(end)>0.05);
        if running==0, messg('Stopping (dchi^2<5 percent)'); end
    else
        messg('Going back to old model and stopping');
        running=0;
        Mod=oldMod;
        dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))    
        CHIQ(end)=[];RMS(end)=[];
        dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
    end
    if CHIQ(end)<1, running=0; end
    if running, % broyden update
        % S = S + ((log(Mod.R)-log(oldR)-S*dM)*(dM'))/(dM'*dM); % BROYDEN
        dM=logtrans(Mod.M(:),lbound,ubound)-logtrans(oldMod.M(:),lbound,ubound);
        if (size(S,1)==length(Mod.R))&&(size(S,2)==length(dM)),
            messg('Doing Broyden update of sensitivity matrix!');
            ddr=(logtrans(Mod.R,lbound,ubound)-logtrans(oldMod.R,lbound,ubound)-S*dM(:))/(dM(:)'*dM(:));
            for i=1:size(S,2), S(:,i)=S(:,i)+ddr*dM(i); end
        end
        if (size(S,2)==length(Mod.R))&&(size(S,1)==length(dM)), %trans+sparse
            messg('Doing Broyden update of sensitivity matrix!');
            ddr=(logtrans(Mod.R,lbound,ubound)-logtrans(oldMod.R,lbound,ubound)-(dM(:)'*S)')'/(dM(:)'*dM(:));
            for i=1:size(S,1), 
                fi=find(S(i,:));
                S(i,fi)=S(i,fi)+ddr(fi)*dM(i);
            end
        end
    end
    % only in full inversion, not in single step
    if isfield(INV,'robust')&&(INV.robust>0),
        dR=(log(N.r)-log(Mod.R))./log(1+N.err);
        w=abs(dR)*sum(abs(dR))/sum(dR.^2);
        w(w<1)=1;
        N.err=N.err.*w;
        messg('Iteratively data reweighting (robust inversion)');
    end    
end


% if ~isfield(Mod,'isfor')||isequal(Mod.isfor,0), 
% %    'Forward response not uptodate! Calculate?','Yes','No','No'
%    dc2dinvres('forwarding_Callback',gcbo,[],guidata(gcbo))
% end
% doit=~isequal(sort(size(S)),sort([length(N.k) (length(Mod.z)-1)*(length(Mod.x)-1)]));
% if doit,
%     dc2dinvres('compsens_Callback',gcbo,[],guidata(gcbo))
% end
% axes(handles.malfeld);
% nM=Mod.M;running=1;dM0=0;
% lbound=0;ubound=0;
% if isfield(INV,'lbound'), lbound=INV.lbound; end
% if isfield(INV,'ubound'), ubound=INV.ubound; end
% if lbound>min(min(N.r),min(Mod.R)), lbound=0; end
% if (ubound>0)&&(ubound<max(max(N.r),max(Mod.R))), ubound=0; end
% while running,
%     % Inverse subproblem with log(Data) and log(Model)
%     if lbound>min(min(Mod.R),min(N.r)), lbound=0;display('resetting lower bound'); end
%     if (ubound>0)&&(ubound<max(max(Mod.R),max(N.r))),
%         ubound=0;display('resetting upper bound'); end
%     if INV.lolo,
%         dR=log(N.r(:)-lbound)-log(Mod.R(:)-lbound);
%         if ubound>0, dR=dR-log(ubound-N.r(:))+log(ubound-Mod.R(:)); end
%         dM0=log(Mod.M(:)-lbound);
%         if isfield(Mod,'Mref')&&isequal(size(Mod.M),size(Mod.Mref)),
%             dM0=dM0-log(Mod.Mref(:)-lbound); 
%             if ubound>0, dM0=dM0-log(ubound-Mod.M(:))+log(ubound-Mod.Mref(:)); end
%         end
%     else
%         dR=N.r(:)-Mod.R(:);
%         if isfield(Mod,'Mref')&&isequal(size(Mod.M),size(Mod.Mref)), 
%             dM0=Mod.M(:)-Mod.Mref(:); end
%     end
%     fi=(INV.lam==0);
%     set(gcf,'Pointer','watch');drawnow;
%     [dM,lam]=invers(S,dR,INV,N.err,dM0);
%     if (INV.auto==1)&&(INV.const==1), 
%         INV.auto=2;INV.lam=lam; 
%     end
%     if fi, INV.first=lam; end
%     % Model updating
%     nM=Mod.M;fak=1; % line search factor
%     if isfield(INV,'linesearch')&&(INV.linesearch>0),
%         if INV.lolo, 
%             if ubound>0,
%                 nM(:)=(ubound*(Mod.M(:)-lbound).*exp(dM)+lbound*(ubound-Mod.M(:)))./((Mod.M(:)-lbound).*exp(dM)+ubound-Mod.M(:));
%             else
%                 nM(:)=(Mod.M(:)-lbound).*exp(dM)+lbound;
%             end
%         else
%             nM(:)=Mod.M(:)+dM; 
%         end
%         R1=abs(dcfwd2d(Mod.x,Mod.z,nM,Mod.Lay,N,FOR));
%         [fak,appR]=linesearch(N,Mod.R,R1,INV.lolo);
%         if fak==0.0,
%             messg('Linear line search failed. Trying parabolic.');
%             taug=0.3; % sampling point
%             nM=Mod.M;
%             if ubound>0,
%                 nM(:)=(ubound*(Mod.M(:)-lbound).*exp(dM*taug)+lbound*(ubound-Mod.M(:)))./((Mod.M(:)-lbound).*exp(dM*taug)+ubound-Mod.M(:));
%             else
%                 nM(:)=(Mod.M(:)-lbound).*exp(dM*taug)+lbound;
%             end
%             Ri=abs(dcfwd2d(Mod.x,Mod.z,nM,Mod.Lay,N,FOR));
%             chiq=[CHIQ(end) chi2(N.r,Ri,N.err,1) chi2(N.r,R1,N.err,1)];
%             taus=[0 taug 1];
%             G=[1 0 0;1 taug taug^2;1 1 1];
%             xx=G\chiq';
%             fak=-xx(2)/2/xx(3);
%         end
%         if fak>0.9, fak=1; end
%         if fak<0, fak=0.1; end
%         if fak<1, dM=dM*fak;
%             messg(sprintf('Line search factor = %.2f chi^2=%.1f-->%.1f',...
%                 fak,chi2(N.r,R1,N.err,INV.lolo),chi2(N.r,appR,N.err,INV.lolo))); end
%     else
%         fak=1.001;
%     end
%     if INV.lolo,
%         if ubound>0,
%            nM(:)=(ubound*(Mod.M(:)-lbound).*exp(dM)+lbound*(ubound-Mod.M(:)))./((Mod.M(:)-lbound).*exp(dM)+ubound-Mod.M(:));
%         else
%             nM(:)=(Mod.M(:)-lbound).*exp(dM)+lbound;
%         end
%     else
%         nM(:)=Mod.M(:)+dM;
%     end
%     %dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
%     axes(handles.malfeld);
% %     draw2dmodel(Mod.x,Mod.z,nM,MAL);
%     patch2dmodel(Mod.x,Mod.z,nM,MAL,N);
%     drawnow;
%     % Forward Calculation
%     oldR=Mod.R;
%     set(gcf,'Pointer','watch');drawnow;
%     if fak==1, Mod.R=R1; elseif fak==0, Mod.R=oldR;
%     else Mod.R=abs(dcfwd2d(Mod.x,Mod.z,nM,Mod.Lay,N,FOR)); end
%     set(gcf,'Pointer','arrow');
%     % stop if rms grows or changes less than 5 percent
%     RMS=[RMS rms(N.r,Mod.R,INV.lolo)];
%     CHIQ=[CHIQ chi2(N.r,Mod.R,N.err,INV.lolo)];
%     dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
%     dchiq=CHIQ(end)-CHIQ(end-1);
%     if (dchiq<0)||(length(RMS)<3),
%         Mod.M=nM;
%         running=(abs(dchiq)/CHIQ(end)>0.05);
%         if running==0, messg('Stopping (dchi^2<5 percent)'); end
%     else
%         messg('Going back to old model and stopping');
%         running=0;
%         Mod.R=oldR;
% %         draw2dmodel(Mod.x,Mod.z,Mod.M,MAL);
%         patch2dmodel(Mod.x,Mod.z,Mod.M,MAL,N);
%         CHIQ(end)=[];RMS(end)=[];
%         dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
%     end
%     if CHIQ(end)<1, running=0; end
%     if running,
%         % S = S + ((log(Mod.R)-log(oldR)-S*dM)*(dM'))/(dM'*dM); % BROYDEN
%         if (size(S,1)==length(Mod.R))&&(size(S,2)==length(dM)),
%             ddr=(log(Mod.R)-log(oldR)-S*dM(:))/(dM(:)'*dM(:));
%             for i=1:size(S,2), S(:,i)=S(:,i)+ddr*dM(i); end
%             messg('Doing Broyden update of sensitivity matrix!');
%         end
%         if (size(S,2)==length(Mod.R))&&(size(S,1)==length(dM)), %trans+sparse
%             ddr=(log(Mod.R)-log(oldR)-(dM(:)'*S)')'/(dM(:)'*dM(:));
%             for i=1:size(S,1), 
%                 fi=find(S(i,:));
%                 S(i,fi)=S(i,fi)+ddr(fi)*dM(i);
%             end
%             messg('Doing Broyden update of sensitivity matrix!');
%         end
%     end
%     if isfield(INV,'robust')&&(INV.robust>0),
%         dR=(log(N.r)-log(Mod.R))./log(1+N.err);
%         w=abs(dR)*sum(abs(dR))/sum(dR.^2);
%         w(w<1)=1;
%         N.err=N.err.*w;
%         messg('Iteratively data reweighting (robust inversion)');
%     end
%     if ishandle(2), dc2dinvres('compare_Callback',gcbo,[],guidata(gcbo)); end
%     dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
%     set(gcf,'Pointer','arrow');
% end
% %dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
% set(handles.backtoold,'Enable','Off');    

% --------------------------------------------------------------------
function varargout = invone_Callback(h, eventdata, handles, varargin)
global Mod N FOR MAL RMS S INV CHIQ
if ~isfield(Mod,'isfor')||isequal(Mod.isfor,0), 
%    'Forward response not uptodate! Calculate?','Yes','No','No'
   dc2dinvres('forwarding_Callback',gcbo,[],guidata(gcbo))
end
set(gcf,'Pointer','watch');drawnow;
save('model.mat','Mod');
doit=~isequal(sort(size(S)),sort([length(N.k) (length(Mod.z)-1)*(length(Mod.x)-1)]));
if doit,
    dc2dinvres('compsens_Callback',gcbo,[],guidata(gcbo))
end
% Inverse subproblem with log(Data-lbound) and log(Model-lbound)
dM0=0;
if INV.lolo,
    lbound=0;ubound=0;
    if isfield(INV,'lbound'), lbound=INV.lbound; end
    if (lbound>min(Mod.R))||(lbound>min(N.r)), lbound=0; end
    if min(Mod.M(:))<lbound, lbound=0; end
    if isfield(INV,'ubound'), ubound=INV.ubound; end
    if ubound<max(max(N.r),max(Mod.R)), ubound=0; end
    if max(Mod.M(:))>ubound, ubound=0; end
    if isfield(Mod,'Mref')&&~isempty(Mod.Mref),
        if isfield(INV,'lbound')&&(lbound>INV.lbound), lbound=0; end %???
        dM0=log(Mod.M(:)-lbound);
        if isfield(Mod,'Mref')&&isequal(size(Mod.M),size(Mod.Mref)),
            dM0=dM0-log(Mod.Mref(:)-lbound); 
            if ubound>0, dM0=dM0-log(ubound-Mod.M(:))+log(ubound-Mod.Mref(:)); end
        end
    end
    dR=log(N.r(:)-lbound)-log(Mod.R(:)-lbound);
    if ubound>0, dR=dR-log(ubound-N.r(:))+log(ubound-Mod.R(:)); end
else
    dR=N.r(:)-Mod.R(:);
    if isfield(Mod,'Mref')&&isequal(size(Mod.M),size(Mod.Mref)), 
        dM0=Mod.M(:)-Mod.Mref(:); end
end
fi=(INV.lam==0);
[dM,lam]=invers(S,dR,INV,N.err,dM0);
if (INV.auto==1)&&(INV.const==1), 
    INV.auto=2;INV.lam=lam; 
end
if fi, INV.first=lam; end
% lb=0;ub=Inf;
lb=0.1;ub=10000;
if min(Mod.M(:))<lb, lb=0; end
if max(Mod.M(:))>ub, ub=Inf; end
nM=Mod.M;fak=1; % line search factor
if isfield(INV,'linesearch')&&(INV.linesearch>0),
    nM05=nM;
    if INV.lolo,
        nM05(:)=(Mod.M(:)-lbound).*exp(dM*0.5)+lbound;
        nM(:)=(Mod.M(:)-lbound).*exp(dM)+lbound;
        if ubound>0,
           nM(:)=(ubound*(Mod.M(:)-lbound).*exp(dM)+lbound*(ubound-Mod.M(:)))./((Mod.M(:)-lbound).*exp(dM)+ubound-Mod.M(:));
           nM05(:)=(ubound*(Mod.M(:)-lbound).*exp(dM*0.5)+lbound*(ubound-Mod.M(:)))./((Mod.M(:)-lbound).*exp(dM*0.5)+ubound-Mod.M(:));
        end
    else
        nM05(:)=Mod.M(:)+dM*0.5;
        nM(:)=Mod.M(:)+dM;
    end
    if INV.linesearch==2,
        Mod.R=abs(dcfwd2d(Mod.x,Mod.z,nM,Mod.Lay,N,FOR));
        R05=abs(dcfwd2d(Mod.x,Mod.z,nM05,Mod.Lay,N,FOR));
        if INV.lolo,
            dR1=log(N.r(:)-lbound)-log(Mod.R(:)-lbound);
            dR05=log(N.r(:)-lbound)-log(R05(:)-lbound);
        else
            dR1=N.r-Mod.R;dR05=N.r-R05;
        end
        Gp=[0 0 1;0.25 0.5 1;1 1 1];GpI=inv(Gp'*Gp);
        y=[norm(dR./N.err);norm(dR05./N.err);norm(dR1./N.err)];
    elseif INV.linesearch==3, % inexact, eigentlich wenig sinnvoll
        Gp=[0 0 1;0.25 0.5 1;1 1 1;4 2 1];GpI=inv(Gp'*Gp);
        dR1=S*dM-dR;dR05=0.5*S*dM-dR;dR2=S*dM*2-dR;
        y=[norm(dR./N.err);norm(dR05./N.err);norm(dR1./N.err);norm(dR2./N.err)];
    else
        R1=abs(dcfwd2d(Mod.x,Mod.z,nM,Mod.Lay,N,FOR));
        [fak,appR]=linesearch(N,Mod.R,R1,INV.lolo);
        if fak==0.0,
            messg('Linear line search failed. Trying parabolic.');
            taug=0.3; % sampling point
            nM=Mod.M;
            if ubound>0,
                nM(:)=(ubound*(Mod.M(:)-lbound).*exp(dM*taug)+lbound*(ubound-Mod.M(:)))./((Mod.M(:)-lbound).*exp(dM*taug)+ubound-Mod.M(:));
            else
                nM(:)=(Mod.M(:)-lbound).*exp(dM*taug)+lbound;
            end
            Ri=abs(dcfwd2d(Mod.x,Mod.z,nM,Mod.Lay,N,FOR));
            chiq=[CHIQ(end) chi2(N.r,Ri,N.err,1) chi2(N.r,R1,N.err,1)];
            taus=[0 taug 1];
            G=[1 0 0;1 taug taug^2;1 1 1];
            xx=G\chiq';
            fak=-xx(2)/2/xx(3);
        end
        if fak>0.9, fak=1; end
        if fak<0, fak=0.1; end
        if fak<1, messg(sprintf('Line search factor = %.2f chi^2=%.1f-->%.1f',...
                fak,chi2(N.r,R1,N.err,INV.lolo),chi2(N.r,appR,N.err,INV.lolo))); end
    end
    if INV.linesearch>1,
        abc=GpI*(Gp'*y);
        xx=0:0.05:2;yy=xx.^2*abc(1)+abc(2)*xx+abc(3);
        figure(1);plot(xx,yy,'-',Gp(:,2),y,'+');
        fak=-abc(2)/abc(1)/2;
        messg(sprintf('Line search factor = %.2f',fak));
    end
end % linesearch specified
running=1;fak0=fak;
while running,
    if INV.lolo,
        if ubound>0,
           nM(:)=(ubound*(Mod.M(:)-lbound).*exp(dM*fak)+lbound*(ubound-Mod.M(:)))./((Mod.M(:)-lbound).*exp(dM*fak)+ubound-Mod.M(:));
        else
            nM(:)=(Mod.M(:)-lbound).*exp(dM*fak)+lbound;
        end
    else
        nM(:)=Mod.M(:)+dM*fak;
    end
    running=0;%(min(nM(:))<lb)|(max(nM(:))>ub);%disabled trusted region
    if running, fak=fak*0.9; end
end
if fak<fak0, messg(sprintf('Trusted region factor %.2f',fak)); end
Mod.M=nM;
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
% Forward Calculation
if isfield(INV,'linesearch')&&(INV.linesearch>0)&&(fak==1), 
    Mod.R=R1; 
else 
    Mod.R=abs(dcfwd2d(Mod.x,Mod.z,Mod.M,Mod.Lay,N,FOR)); 
end
RMS=[RMS rms(N.r,Mod.R,INV.lolo)];
CHIQ=[CHIQ chi2(N.r,Mod.R,N.err,INV.lolo)];
dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
if ishandle(2), dc2dinvres('compare_Callback',gcbo,[],guidata(gcbo)); end
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
set(gcf,'Pointer','arrow');
set(handles.backtoold,'Enable','On');    

% -------------------------- SHOW ------------------------------------
% --------------------------------------------------------------------
function varargout = show_Callback(h, eventdata, handles, varargin)
global Mod N S
if isfield(N,'r'), lr=length(N.r); else lr=0; end
ss=size(S);if ss(2)==lr, ss=fliplr(ss); end
stat='Off';if isfield(Mod,'M')&&(size(S,2)==numel(Mod.M)), stat='On'; end
set(handles.showcov,'Enable',stat);
set(handles.showsens,'Enable',stat);
if (lr>0)&&isequal(sort(size(S)),sort([lr numel(Mod.M)])),
    set(handles.showsens,'Enable','On'); else
    set(handles.showsens,'Enable','Off'); end
if (lr>0)&&isfield(N,'u')&&(length(N.u)==lr)
    set(handles.showvoltage,'Enable','On'); else
    set(handles.showvoltage,'Enable','Off'); end
if (lr>0)&&isfield(N,'i')&&(length(N.i)==lr)
    set(handles.showcurrent,'Enable','On'); else
    set(handles.showcurrent,'Enable','Off'); end
if (lr>0)&&isfield(N,'ip')&&(length(N.ip)==lr)
    set(handles.showip,'Enable','On'); else
    set(handles.showip,'Enable','Off'); end
if (lr>0)&&isfield(N,'ip')&&isfield(Mod,'IP')&&(length(N.ip)==lr)&&(ss(2)>0)&&(ss(2)==numel(Mod.IP)),
    set(handles.compip,'Enable','On'); else
    set(handles.compip,'Enable','Off'); end
if (lr>0)&&isfield(Mod,'IP')&&isequal(size(Mod.IP),size(Mod.M))&&(max(size(Mod.IP))>0),
    set(handles.showcharge,'Enable','On'); else
    set(handles.showcharge,'Enable','Off'); end
set(handles.showmod,'Enable',sonoff(isfield(Mod,'M')&&~isempty(Mod.M)));
ssens=sonoff(isfield(Mod,'M')&&~isempty(Mod.M)&&(ss(2)==numel(Mod.M)));
set(handles.showsens,'Enable',ssens);
set(handles.showcov,'Enable',ssens);
isdata=isfield(N,'a')&isfield(N,'r');
set(handles.showdata,'Enable',sonoff(isdata));
isres=(isdata&&(length(N.r)==length(Mod.R)));
set(handles.showfor,'Enable',sonoff(isres));
set(handles.showmisfit,'Enable',sonoff(isres));
set(handles.compare,'Enable',sonoff(isres));
iserr=isdata&&isfield(N,'err')&&(length(N.err)==length(N.r));
set(handles.showerror,'Enable',sonoff(iserr));
set(handles.misnorm,'Enable',sonoff(isres&iserr));

% --------------------------------------------------------------------
function varargout = showsens_Callback(h, eventdata, handles, varargin)
global S malstat %Mod.x Mod.z MAL
set(handles.preview,'Visible','Off');
set(handles.slider,'Visible','On');
set(handles.editslider,'Visible','On');
set(handles.dispslider,'Visible','On');
set(handles.dispslider,'String','Show Sensitivity for datum');

ma=size(S,1);
set(handles.slider,'sliderstep',[1 1]*(1/(ma-1)),'max',ma,'min',1,'Value',1);
set(handles.editslider,'String','1');
malstat=9;
dc2dinvres('slider_Callback',gcbo,[],guidata(gcbo))
set(handles.widerstand,'Visible','Off');
set(handles.widerstand,'String','0');
set(handles.ipmodelinput,'Visible','Off');
set(handles.ipmodelinput,'String','0');
%set(handles.helptext,'Visible','Off');
%'String','Click to see/delete data');

% --------------------------------------------------------------------
function varargout = showcov_Callback(h, eventdata, handles, varargin)
global Mod S N malstat
axes(handles.malfeld);
ma=struct('clog',1,'cmap',0);
% Mod.Cov=zeros(size(Mod.M));
% if size(S,2)Mod.Cov(:)=sum(abs(S));
% draw2dmodel(Mod.x,Mod.z,log10(Mod.Cov),ma);
patch2dmodel(Mod.x,Mod.z,Mod.Cov,ma,N);
malstat=4;
set(handles.preview,'Visible','Off');
set(handles.slider,'Visible','Off');
set(handles.editslider,'Visible','Off');
set(handles.dispslider,'Visible','Off');
set(handles.widerstand,'Visible','Off');
set(handles.ipmodelinput,'Visible','Off');
set(handles.helptext,'Visible','On','String','Coverage');

% --------------------------------------------------------------------
function varargout = showdata_Callback(h, eventdata, handles, varargin)
global N malstat MAL
axes(handles.malfeld);
mal=MAL;mal.canot='Ohm*m';
showdata2d(N,N.r,mal);
malstat=1;
im=get(handles.malfeld,'Children');
set(im,'ButtonDownFcn','dc2dinvres(''malfeld_ButtonDownFcn'',gcbo,[],guidata(gcbo))')
set(handles.preview,'Visible','Off');
set(handles.slider,'Visible','Off');
set(handles.editslider,'Visible','Off');
set(handles.dispslider,'Visible','Off');
set(handles.helptext,'String','RhoA in OhmM - Click to see/delete data');
set(handles.helptext,'Visible','On');
set(handles.widerstand,'Visible','Off');
set(handles.ipmodelinput,'Visible','Off');


% --------------------------------------------------------------------
function showerror_Callback(hObject, eventdata, handles)
global N malstat MAL
if ~isfield(N,'err'), return; end
mal=MAL;mal.cauto=1;mal.canot='err/%';mal.clog=1;
axes(handles.malfeld);
showdata2d(N,N.err*100,mal);
malstat=-1;
im=get(handles.malfeld,'Children');
set(im,'ButtonDownFcn','dc2dinvres(''malfeld_ButtonDownFcn'',gcbo,[],guidata(gcbo))')
set(handles.preview,'Visible','Off');
set(handles.slider,'Visible','Off');
set(handles.editslider,'Visible','Off');
set(handles.dispslider,'Visible','Off');
set(handles.helptext,'String','Error in % Click to see/delete data');
set(handles.widerstand,'Visible','Off');
set(handles.ipmodelinput,'Visible','Off');

% --------------------------------------------------------------------
function misnorm_Callback(hObject, eventdata, handles)
global N malstat MAL Mod
mal=MAL;mal.cauto=1;mal.clog=0;mal.cmap=2;
% mal.canot='misfit/error';
axes(handles.malfeld);
dR=1-Mod.R./N.r;
if isfield(N,'err'), showdata2d(N,dR./N.err,mal); end
malstat=-1;
%im=get(handles.malfeld,'Children');
%set(im,'ButtonDownFcn','dc2dinvres(''malfeld_ButtonDownFcn'',gcbo,[],guidata(gcbo))')
set(handles.preview,'Visible','Off');
set(handles.slider,'Visible','Off');
set(handles.editslider,'Visible','Off');
set(handles.dispslider,'Visible','Off');
set(handles.helptext,'String','Misfit/Error: Click to see/delete data');
set(handles.widerstand,'Visible','Off');
set(handles.ipmodelinput,'Visible','Off');

% --------------------------------------------------------------------
function varargout = showfor_Callback(h, eventdata, handles, varargin)
global N malstat Mod MAL
axes(handles.malfeld);
showdata2d(N,Mod.R,MAL);
malstat=2;
set(handles.preview,'Visible','Off');
set(handles.slider,'Visible','Off');
set(handles.editslider,'Visible','Off');
set(handles.dispslider,'Visible','Off');
set(handles.helptext,'String','');
set(handles.widerstand,'Visible','Off');
set(handles.ipmodelinput,'Visible','Off');

% --------------------------------------------------------------------
function varargout = showmisfit_Callback(h, eventdata, handles, varargin)
global N Mod malstat MAL %INV
mal=struct('cauto',1,'cmap',2,'clog',0);
if isfield(MAL,'xdir'), mal.xdir=MAL.xdir; end
axes(handles.malfeld);
%if INV.lolo,
%    dR=(log(N.r(:))./log(Mod.R(:))-1)*100;
%else
dR=(1-Mod.R(:)./N.r(:))*100;
%end
showdata2d(N,dR,mal);
malstat=3;
set(handles.preview,'Visible','Off');
set(handles.slider,'Visible','Off');
set(handles.editslider,'Visible','Off');
set(handles.dispslider,'Visible','Off');
set(handles.helptext,'Visible','On',...
    'String','Data misfit in %');
set(handles.widerstand,'Visible','Off');
set(handles.ipmodelinput,'Visible','Off');

% --------------------------------------------------------------------
function varargout = compare_Callback(h, eventdata, handles, varargin)
global N Mod MAL
%dR=(log(Mod.R(:))-log(N.r(:)))./log(N.r(:))*100;
dR=(1-Mod.R(:)./N.r(:))*100;
figure(2);
set(2,'MenuBar','none','NumberTitle','off','Name','Data,Model response & misfit');
iconify(2);
mal=MAL;mal.cauto=0;
mal.cmax=max([N.r(:);Mod.R(:)]);
mal.cmin=min([N.r(:);Mod.R(:)]);
subplot(3,1,1);
showdata2d(N,[],mal);
xl=get(gca,'XLim');xl=xl(1)+0.05*diff(xl);
yl=get(gca,'YLim');yl=yl(1)+0.9*diff(yl);
text(xl,yl,'Measured Data in \Omega m');
subplot(3,1,2);
showdata2d(N,Mod.R,mal);
text(xl,yl,'Calculated Data in \Omega m');
subplot(3,1,3);
mal.log=0;mal.cauto=0;
mm=max(abs(dR));mal.cmin=-mm;mal.cmax=mm;
showdata2d(N,dR,mal);
%colormap(jet);
text(xl,yl,'Data Misfit in %');

% --------------------------------------------------------------------
function varargout = showdoi_Callback(h, eventdata, handles, varargin)
global Mod N S INV malstat
axes(handles.malfeld);
% draw2dmodel(Mod.x,Mod.z,getdoi(N.r,S,Mod.M(1,1),1.05,INV));
patch2dmodel(Mod.x,Mod.z,getdoi(N.r,S,Mod.M(1,1),1.05,INV),[],N);
malstat=5;
set(handles.preview,'Visible','Off');
set(handles.slider,'Visible','Off');
set(handles.editslider,'Visible','Off');
set(handles.dispslider,'Visible','Off');
set(handles.helptext,'String','');
set(handles.widerstand,'Visible','Off');
set(handles.ipmodelinput,'Visible','Off');

% --------------------------------------------------------------------
function varargout = showmod_Callback(h, eventdata, handles, varargin)
global Mod MAL malstat N
if ishandle(2), dc2dinvres('compare_Callback',gcbo,[],guidata(gcbo)); end
axes(handles.malfeld);
mal=MAL;mal.canot='Ohm*m';
if isfield(MAL,'style')&&(MAL.style==3),
    set(figure(8),'MenuBar','none','NumberTitle','off','Name','Reciprocity crossplot');
    plotlongmod(Mod,mal);
    return;
end
if isfield(MAL,'alpha')&&(MAL.alpha)&&isfield(Mod,'Cov'), 
%         draw2dmodel(Mod.x,Mod.z,Mod.M,MAL,Mod.Cov);
        patch2dmodel(Mod.x,Mod.z,Mod.M,mal,N,Mod.Cov);
else
%         draw2dmodel(Mod.x,Mod.z,Mod.M,MAL);
        patch2dmodel(Mod.x,Mod.z,Mod.M,mal,N);
end
global FIX XX ZZ
xz=zeros(size(Mod.x));zz=Mod.z;
if isfield(N,'topo'),
    xz=interp1(N.topo(:,1),N.topo(:,2),Mod.x,'linear','extrap');
    zz=-Mod.z;
end

hold on
if isequal(size(FIX),size(Mod.M)),
    [i,j]=find(FIX==-1);
    plot((Mod.x(i)+Mod.x(i+1))/2,(zz(j)+zz(j+1))/2+xz(i),'wx');
    [i,j]=find(FIX>0);
    plot((Mod.x(i)+Mod.x(i+1))/2,(zz(j)+zz(j+1))/2+xz(i),'wo');
end
if isequal(size(XX)+[1 0],size(Mod.M)),
    [i,j]=find(XX>0);
    plot(Mod.x(i+1),(zz(j)+zz(j+1))/2+xz(i),'w+');
end
if isequal(size(ZZ)+[0 1],size(Mod.M)),
    [i,j]=find(ZZ>0);
    plot((Mod.x(i)+Mod.x(i+1))/2,zz(j+1)+xz(i),'w+');    
end
hold off
malstat=0;
im=get(handles.malfeld,'Children');
set(im,'ButtonDownFcn','dc2dinvres(''malfeld_Model'',gcbo,[],guidata(gcbo))')
%set(im,'ButtonUpFcn','dc2dinvres(''malfeldup_Model'',gcbo,[],guidata(gcbo))')
set(handles.preview,'Visible','Off');
set(handles.slider,'Visible','Off');
set(handles.editslider,'Visible','Off');
set(handles.dispslider,'Visible','Off');
set(handles.helptext,'Visible','On',...
    'String',{'Resistivity in Ohmm, click to see','type resistivity to change cells','F/U for fixing/unfixing resistivities'});
set(handles.widerstand,'Visible','On');
if isfield(N,'ip')&&(length(N.ip)==length(N.r)), set(handles.ipmodelinput,'Visible','On'); end
% set(handles.widerstand,'String','');
set(gcf,'Pointer','arrow');
drawnow;

% --------------------------------------------------------------------
function varargout = dc2dinvres_ResizeFcn(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = dc2dinvres_KeyPressFcn(h, eventdata, handles, varargin)
aa=get(gcf,'CurrentCharacter');
set(gcf,'Pointer','arrow');
switch aa, %K N Y Z
    case 'A'
        dc2dinvres('addlb_Callback',gcbo,[],guidata(gcbo))
    case 'B'
        dc2dinvres('blockify_Callback',gcbo,[],guidata(gcbo))
    case 'C'
        dc2dinvres('compare_Callback',gcbo,[],guidata(gcbo))
    case 'D'
        dc2dinvres('showdata_Callback',gcbo,[],guidata(gcbo))
    case 'E'
        dc2dinvres('fileexport_Callback',gcbo,[],guidata(gcbo))
    case 'F'
        dc2dinvres('forwarding_Callback',gcbo,[],guidata(gcbo))
    case 'G'
        dc2dinvres('optgraf_Callback',gcbo,[],guidata(gcbo))
    case 'H'
        dc2dinvres('starthom_Callback',gcbo,[],guidata(gcbo))
    case 'I'
        dc2dinvres('invauto_Callback',gcbo,[],guidata(gcbo))
    case 'J'
        dc2dinvres('compsens_Callback',gcbo,[],guidata(gcbo))
    case '1'
        dc2dinvres('invone_Callback',gcbo,[],guidata(gcbo))
    case 'L'
        dc2dinvres('fileload_Callback',gcbo,[],guidata(gcbo))
    case 'M'
        dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
    case 'O'
        dc2dinvres('fileread_Callback',gcbo,[],guidata(gcbo))
    case 'P'
        dc2dinvres('showcharge_Callback',gcbo,[],guidata(gcbo))
    case 'Q'
        dc2dinvres('readuline_Callback',gcbo,[],guidata(gcbo))
    case 'R'
        dc2dinvres('showfor_Callback',gcbo,[],guidata(gcbo))
    case 'S'
        dc2dinvres('filesave_Callback',gcbo,[],guidata(gcbo))
    case 'T'
        dc2dinvres('timelapse_Callback',gcbo,[],guidata(gcbo))
    case 'U'
        dc2dinvres('readline_Callback',gcbo,[],guidata(gcbo))
    case 'V'
        dc2dinvres('vtkexport_Callback',gcbo,[],guidata(gcbo))
    case 'W'
        dc2dinvres('filewrite_Callback',gcbo,[],guidata(gcbo))
    case 'X'
        dc2dinvres('fileexit_Callback',gcbo,[],guidata(gcbo))
    case 'backspace'
        dc2dinvres('backtoold_Callback',gcbo,[],guidata(gcbo))
end
% stillfree= K N Y Z

% --------------------------------------------------------------------
function varargout = forwarding_Callback(h, eventdata, handles, varargin)
global Mod N FOR INV RMS CHIQ S
set(gcf,'Pointer','watch');drawnow;
if max(N.elec(:,2))>0, %problem with 1d and subsurface electrodes
    Mod.M(end,end)=Mod.M(end,end)*1.0001;
    Mod.R=dcfwd2d(Mod.x,Mod.z,Mod.M,Mod.Lay,N,FOR);
else   
    islay=1;lay=0;
    for i=1:size(Mod.M,2),
        lay(i)=max(Mod.M(:,i));
        islay=islay&(lay(i)==min(Mod.M(:,i)));
    end
    if islay, % layering detected
        if length(unique(lay))==1, %constant model
            Mod.R=ones(size(N.r))*lay(1);
        else
            lay(end+1)=lay(end); % last layer
            messg('Calculating forward response of 1D model...');
            Mod.R=fwd2d1d(N,lay,Mod.z);
        end
    else
        Mod.R=dcfwd2d(Mod.x,Mod.z,Mod.M,Mod.Lay,N,FOR);
    end
end
Mod.isfor=1;
if isfield(Mod,'IP')&&isequal(sort(size(S)),sort([numel(Mod.IP) length(N.a)])),
    if size(S,2)==numel(Mod.IP), Mod.ipfor=S*Mod.IP(:); end
    if size(S,1)==numel(Mod.IP), Mod.ipfor=(Mod.IP(:)'*S)'; end    
end
RMS=[RMS rms(N.r,Mod.R,INV.lolo)];
CHIQ=[CHIQ chi2(N.r,Mod.R,N.err,INV.lolo)];
dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
set(gcf,'Pointer','arrow');
if ishandle(2), dc2dinvres('compare_Callback',gcbo,[],guidata(gcbo))
else dc2dinvres('showfor_Callback',gcbo,[],guidata(gcbo)); end
% if isfield(N,'ip')&&(length(N.ip)==length(N.r)),
%     global S IP
% end

% --------------------------------------------------------------------
function varargout = modpar_Callback(h, eventdata, handles, varargin)
global Mod RMS FOR INV N CHIQ
s_mod;
messg(strcat(sprintf('Model(%dx%d=%d cells): Mod.x=%g..%g',size(Mod.M,1),size(Mod.M,2),numel(Mod.M),min(Mod.x),max(Mod.x)),...
    '  Mod.z= ',sprintf('%g ',Mod.z)));  
dc2dinvres('loadsens_Callback',gcbo,[],guidata(gcbo))
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
if find(diff(Mod.Lay)), % apparently a layered model
    set(gcf,'Pointer','watch');drawnow;
    if max(N.elec(:,2))>0,
        Mod.R=dcfwd2d(Mod.x,Mod.z,Mod.M,Mod.Lay,N,FOR);
    else
        messg('Calculating forward response of 1D model...');
        Mod.R=fwd2d1d(N,Mod.Lay,Mod.z);
    end
    Mod.Mref=Mod.M;
    set(gcf,'Pointer','arrow');
else
    Mod.R=ones(size(N.r))*Mod.Lay(1);
end
if ~isequal(size(Mod.M),size(Mod.Mref)), Mod.Mref=Mod.M; end
RMS=rms(N.r,Mod.R,INV.lolo);
CHIQ=chi2(N.r,Mod.R,N.err,INV.lolo);
dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
if ishandle(2), dc2dinvres('compare_Callback',gcbo,[],guidata(gcbo)); end

% --------------------------------------------------------------------
function varargout = modeledit_Callback(h, eventdata, handles, varargin)
disp('Not yet implemented!');


% SVD ANALYSIS MENU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
function varargout = svdana_Callback(h, eventdata, handles, varargin)
global S RM Mod %RD
onoff='Off';
nm=size(S,2);
if (nm>0)&&(nm==size(RM,1)), onoff='On'; end
set(handles.modcov,'Enable',onoff);
set(handles.modres,'Enable',onoff);
set(handles.datres,'Enable',onoff);
set(handles.datcov,'Enable',onoff);
set(handles.resrad,'Enable',onoff);
set(handles.filfak,'Enable',onoff);
set(handles.sv,'Enable',onoff);
if isfield(Mod,'M')&&(nm==numel(Mod.M)), set(handles.modcov,'Enable','on'); end

% --------------------------------------------------------------------
function varargout = sv_Callback(h, eventdata, handles, varargin)
global s VD %Mod N VM S
if isempty(VD),
    dc2dinvres('compsvd_Callback',gcbo,[],guidata(gcbo))
end
figure(1);
set(1,'MenuBar','none','NumberTitle','off','Name','Singular values');
semilogy(s);

% --------------------------------------------------------------------
function varargout = modcov_Callback(h, eventdata, handles, varargin)
global Mod N s VD VM S MAL RM malstat
if isempty(VD),
%     dc2dinvres('compsvd_Callback',gcbo,[],guidata(gcbo))
end
set(handles.helptext,'String','');
set(handles.slider,'Visible','On');
set(handles.editslider,'Visible','On');
set(handles.dispslider,'Visible','On');
set(handles.dispslider,'String','Covariance for model cell');
ma=size(S,2);
set(handles.slider,'sliderstep',[1 1]*(1/(ma-1)),'max',ma,'min',1,'Value',round((length(Mod.x)-1)*3.5));
set(handles.editslider,'String','1');
malstat=11;
dc2dinvres('slider_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function varargout = datcov_Callback(h, eventdata, handles, varargin)
global Mod N malstat s VD VM S RD
if isempty(VD),
    dc2dinvres('compsvd_Callback',gcbo,[],guidata(gcbo))
end
set(handles.helptext,'String','');
set(handles.slider,'Visible','On');
set(handles.editslider,'Visible','On');
set(handles.dispslider,'Visible','On');
set(handles.dispslider,'String','Covariance for datum');
ma=size(RD,2);
set(handles.slider,'sliderstep',[1 1]*(1/(ma-1)),'max',ma,'min',1,'Value',1);
set(handles.editslider,'String','1');
axes(handles.malfeld);
malstat=10;
dc2dinvres('slider_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function varargout = modres_Callback(h, eventdata, handles, varargin)
global Mod N MAL malstat s VD VM S RM RD
if isempty(VD),
    dc2dinvres('compsvd_Callback',gcbo,[],guidata(gcbo))
end
Rm=Mod.M;
Rm(:)=diag(RM);
messg(sprintf('Model information %.1f, mean %d%%',sum(Rm(:)),round(mean(Rm(:))*100)));
mal=MAL;
mal.cauto=0;mal.log=0;mal.cmin=0;mal.cmax=1;
axes(handles.malfeld);
if 1, %normal
%     draw2dmodel(Mod.x,Mod.z,Rm,mal);
    patch2dmodel(Mod.x,Mod.z,Rm,mal,N);
else %with distortion flag
    mRm=Mod.M;
    mRm(:)=max(RM);
    mRm=Rm./mRm;
    mRm(mRm<0)=0;
    mal.noascal=1;
%     draw2dmodel(Mod.x,Mod.z,Rm,mal,mRm);
    patch2dmodel(Mod.x,Mod.z,Rm,mal,mRm);
end
malstat=6;
set(handles.preview,'Visible','Off');
set(handles.slider,'Visible','Off');
set(handles.editslider,'Visible','Off');
set(handles.dispslider,'Visible','Off');
set(handles.helptext,'String','');

% --------------------------------------------------------------------
function varargout = datres_Callback(h, eventdata, handles, varargin)
global Mod N MAL s VD VM S malstat RM RD
if isempty(VD),
    dc2dinvres('compsvd_Callback',gcbo,[],guidata(gcbo))
end
axes(handles.malfeld);
mal=MAL;
mal.cauto=0;mal.log=0;mal.cmin=0;mal.cmax=1;
Rd=diag(RD);
messg(sprintf('Information content %.1f, effectivity %d%%',sum(Rd),round(mean(Rd)*100)));
showdata2d(N,Rd,mal);
malstat=7;
set(handles.preview,'Visible','Off');
set(handles.slider,'Visible','Off');
set(handles.editslider,'Visible','Off');
set(handles.dispslider,'Visible','Off');
set(handles.helptext,'String','');

% --------------------------------------------------------------------
function varargout = resrad_Callback(h, eventdata, handles, varargin)
global Mod N MAL s VD VM S malstat RM RD
if isempty(VD),
    dc2dinvres('compsvd_Callback',gcbo,[],guidata(gcbo))
end
Rm=abs(diag(RM));
axes(handles.malfeld);
[DX,DZ]=ndgrid(diff(Mod.x),diff(Mod.z));
Fl=DX.*DZ;
Resrad=zeros(size(Mod.M));
Resrad(:)=sqrt(Fl(:)./Rm(:)/pi); % after Friedel (boxcar)
%Resrad(:)=sqrt(Fl(:))/2./erfinv(Rm/2); % equivalent Gaussian stdev
mal=MAL;mal.cauto=1;mal.log=1;
% draw2dmodel(Mod.x,Mod.z,Resrad,mal);
patch2dmodel(Mod.x,Mod.z,Resrad,mal,N);
malstat=8;
set(handles.preview,'Visible','Off');
set(handles.slider,'Visible','Off');
set(handles.editslider,'Visible','Off');
set(handles.dispslider,'Visible','Off');
set(handles.helptext,'String','');

% --------------------------------------------------------------------
function varargout = compsens_Callback(h, eventdata, handles, varargin)
global S Mod N datfile INV
set(gcf,'Pointer','watch');
S=calcsens2d(Mod.x,Mod.z,N,INV);
Mod.Cov=zeros(size(Mod.M));
nm=numel(Mod.M);
if size(S,2)==nm, 
%     Mod.Cov(:)=sum(abs(S)); 
    for i=1:nm, Mod.Cov(i)=sum(abs(S(:,i))); end
end
if size(S,1)==nm, 
%     Mod.Cov(:)=sum(abs(S),2); 
    for i=1:nm, Mod.Cov(i)=sum(abs(S(i,:))); end
end
set(gcf,'Pointer','arrow');
if ~isempty(datfile),
    [fpath,name,ext]=fileparts(datfile);
    matfile=strrep(datfile,ext,'-sens.mat');
    xsave=Mod.x;zsave=Mod.z;
    if ispc&&(exist(matfile,'file')),
        dos(['attrib -H +A "' matfile '"']);
%         fprintf('Setting attrib...');
    end
    save(matfile,'S','xsave','zsave');
    %     fileattrib(matfile,'+h -a');
    if ispc, dos(['attrib +H -A "' matfile '"']); end
end
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
previewon(handles);

% --------------------------------------------------------------------
function varargout = compsvd_Callback(h, eventdata, handles, varargin)
global VD VM s RM RD N datfile Mod INV S
si=[length(N.a) (length(Mod.x)-1)*(length(Mod.z)-1)];
if ~isequal(size(S),si),
    error('Size of S mismatches data/model length!');
end
if INV.method==2, % SIRT
    cov=sum(abs(S));G=spdiags(1./cov(:),0,si(2),si(2));
    voc=sum(abs(S),2);E=spdiags(voc(:),0,si(1),si(1));
    % eigentlich 1./voc(:), aber es funzt!!!
    RM=(S*G)'*E*S;
    RD=S*G*S'*E;
%     aa=sum(RM.^2);
%     RM=RM*spdiags(1./sqrt(aa(:)),0,length(aa),length(aa));
%     bb=sum(RD.^2);
%     RD=RD*spdiags(1./sqrt(bb(:)),0,length(bb),length(bb));
%     error('No resolution available for SIRT!');%doch!
    return; 
end
[fpath,name,ext]=fileparts(datfile);
matfile=strrep(datfile,ext,'-sens.mat');

% if exist(matfile)==2,
%     load(matfile);
%     if (~isequal(Mod.x,xsave))|(~isequal(Mod.z,zsave)),
%         global S
%         dc2dinvres('compsens_Callback',gcbo,[],guidata(gcbo));
%     end
% else
%     global S 
%     dc2dinvres('compsens_Callback',gcbo,[],guidata(gcbo))
% end
t0=clock;
data=length(N.r);
if ~isfield(N,'err'),
    daterr;
    messg(sprintf('Error min=%.1f%% max=%.1f%% mean=%.1f%%',...
        min(N.err)*100,max(N.err)*100,mean(N.err)*100));
end
D=spdiags(1./log(1+N.err),0,data,data);
set(gcf,'Pointer','watch');
nsvd=(INV.method==1)|(INV.weight==0); %normal svd (tsvd & marquardt)
if nsvd,
    messg('Computing singular value decomposition...');
    [VD,s,VM]=svd(D*S);  % Model & Data Vectors
    s=diag(s);
else
    messg('Computing generalized singular value decomposition...');
    if INV.weight==1, C=smooth2d1st(Mod.x,Mod.z); end
    if INV.weight==2, C1=smooth2d2nd(Mod.x,Mod.z);C=C1'*C1; end
    if INV.weight==3, lc=numel(Mod.Cov);C=spdiags(1./Mod.Cov(:),0,lc,lc); end
    [VD,sm,VM]=gensvd(D*S,C);s=sm(:,1)./sm(:,2);
end
set(gcf,'Pointer','arrow');
lambda=0;
if isfield(INV,'first')&&(INV.first>0),
    lambda=INV.first;
elseif isfield(INV,'lam')&&(INV.lam>0),
    lambda=INV.lam;
else
    dR=log(N.r)-log(Mod.Lay(1));
    [dM,lambda]=invshifted(D*S,D*dR,100,0.001,0.5); %???
end
lambda=sqrt(lambda);
messg(sprintf('ready(%.2fs) min(sv)= %g, max(sv)= %g',...
    etime(clock,t0),min(s),max(s)));
% Regularization
r=max(find(s));
if INV.method==1,%tsvd
    if isfield(INV,'first'), r=INV.first; end
    f=ones(r,1);
else
    f=s(1:r).^2./(s(1:r).^2+lambda^2);
end
RD=VD(:,1:r)*diag(f)*VD(:,1:r)';
if nsvd,
    RM=VM(:,1:r)*diag(f)*VM(:,1:r)';
else
    ff=[f;zeros(size(VM,1)-length(f),1)];
    RM=VM*diag(ff)*inv(VM);
end

% --------------------------------------------------------------------
function filfak_Callback(hObject, eventdata, handles)
global s INV
if isempty(s),
    dc2dinvres('compsvd_Callback',gcbo,[],guidata(gcbo));
end
if isfield(INV,'first')&&(INV.first>0), filfak(s,INV.first);
elseif isfield(INV,'lam')&&(INV.lam>0), filfak(s,INV.first); end

% --------------------------------------------------------------------
function varargout = clearsvd_Callback(h, eventdata, handles, varargin)
global VD VM s RM RD malstat
s=[];VD=[];VM=[];RM=[];RD=[];
if ismember(malstat,[6 7 8 10 11]),
    dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo)); end

% --------------------------------------------------------------------
function varargout = clearsens_Callback(h, eventdata, handles, varargin)
global S
S=[];Mod.Cov=[];

% --------------------------------------------------------------------
function varargout = exportmodel_Callback(h, eventdata, handles, varargin)
global Mod datfile S N
[fpath,name,ext]=fileparts(datfile);
if strcmp(ext,''),
    outfile='*.mod';
else
    outfile=strrep(datfile,ext,'.mod');
end
[fname,pname]=uiputfile(outfile,'Export model to ASCII-File');
if fname~=0,
    outfile=fullfile(pname,fname);
    [fpath,name,ee]=fileparts(outfile);
    if strcmp(ee,''), outfile=[outfile '.mod']; end
    messg(['Exporting model to ASCII File : ' outfile]);
    sec=Mod.Cov;isip=0;
    if isfield(Mod,'IP')&&isequal(size(Mod.M),size(Mod.IP)), 
        sec=Mod.IP;isip=1; end
%     modelexport2d(outfile,Mod.M,Mod.x,Mod.z,IP,0,zeros(size(Mod.x)),1);
    if isfield(N,'topo'),
        xz=interp1(N.topo(:,1),N.topo(:,2),Mod.x,'linear','extrap');
    else
        xz=zeros(size(Mod.x));
    end
    %em=questdlg('Saving by midpoint or corner points?','Export Mode','Midpoint','Corner','Midpoint');
    modelexport2d(outfile,Mod.M,Mod.x,Mod.z,sec,0,xz,isip);
    %             modelexport2d(outfile,Mod.M,Mod.x,Mod.z,Mod.Cov,0,xz); % ,isequal(em,'Corner')
    %         else
    %             modelexport2d(outfile,Mod.M,Mod.x,Mod.z,Mod.Cov); % ,isequal(em,'Corner')
    %         end
% end
%     warndlg('Resistivity values are replaced by zero in test version!','Warning');
end % Testversion *0, + Zeile

% ---------------------------MOUSE-BUTTON-----------------------------------
function varargout = malfeld_ButtonDownFcn(h, eventdata, handles, varargin)
global N malstat RMS CHIQ Mod INV S
cp=get(handles.malfeld,'CurrentPoint');
[mids,konfs,ii,kk]=midkonf2d(N);
%[aa,nk]=min(abs(konfs-cp(1,2)));
nk=round(cp(1,2));
[aa,nx]=min(abs(mids-cp(1,1)));
nr=find((kk==nk)&(ii==nx));
if ~isempty(nr),
    if length(nr)>1, 
        %         fprintf('%d \n',length(nr)); 
        nr=nr(end);
    end
    mess=sprintf('Nr %d: rhoa=%.1f',nr,N.r(nr));
    if isfield(N,'ip')&&(~isempty(N.ip))&&(length(N.ip)>=nr), 
        mess=sprintf('%s Phase=%.2fmrad',mess,N.ip(nr)); end  
    if isfield(N,'i')&&(length(N.i)>=nr), 
        mess=sprintf('%s I=%.1fmA',mess,N.i(nr)*1000); end  
    if isfield(N,'u')&&(length(N.u)>=nr), 
        mess=sprintf('%s U=%.1fmV',mess,N.u(nr)*1000); end  
    if isfield(N,'err')&&(~isempty(N.err))&&(length(N.err)>=nr), 
        mess=sprintf('%s Error=%.1f%%',mess,N.err(nr)*100); end  
    mess=sprintf('%s\na=%d (%.1fm)',mess,N.a(nr),N.elec(N.a(nr),1));
    if N.b(nr)>0, mess=sprintf('%s  b=%d (%.1fm)',mess,N.b(nr),N.elec(N.b(nr),1)); end
    mess=sprintf('%s  m=%d (%.1fm)',mess,N.m(nr),N.elec(N.m(nr),1));
    if N.n(nr)>0, mess=sprintf('%s  n=%d (%.1fm)',mess,N.n(nr),N.elec(N.n(nr),1)); end
    if strcmp(questdlg(mess,'Delete Datum?','Yes','No','No'),'Yes'),
        N.a(nr)=[];N.b(nr)=[];N.m(nr)=[];N.n(nr)=[];
        N.r(nr)=[];N.k(nr)=[];
        if isfield(N,'err')&&(length(N.err)>=nr), N.err(nr)=[]; end
        if isfield(N,'ip')&&(length(N.ip)>=nr), N.ip(nr)=[]; end
        if isfield(N,'i')&&(length(N.i)>=nr), N.i(nr)=[]; end
        if isfield(N,'u')&&(length(N.u)>=nr), N.u(nr)=[]; end
        if isfield(N,'konf')&&(length(N.konf)>=nr), N.konf(nr)=[]; end
        if isfield(N,'rho')&&(length(N.rho)>=nr), N.rho(nr)=[]; end        
        Mod.R(nr)=[];
        if (size(S,1)==length(N.r))&&(size(S,1)>=nr), S(nr,:)=[]; end
        if (size(S,2)==length(N.r))&&(size(S,2)>=nr), S(:,nr)=[]; end
        RMS(end)=rms(N.r,Mod.R,INV.lolo);
        CHIQ(end)=chi2(N.r,Mod.R,N.err,INV.lolo);
        messg(sprintf('Deleting Datum Number %d from data set. new Chi^2=%.1f RMS=%.1f',nr,CHIQ(end),RMS(end)));
        if malstat==1, dc2dinvres('showdata_Callback',gcbo,[],guidata(gcbo)); end
        if malstat==-1, dc2dinvres('showerror_Callback',gcbo,[],guidata(gcbo)); end
        if malstat==12, dc2dinvres('showvoltage_Callback',gcbo,[],guidata(gcbo)); end
        if malstat==13, dc2dinvres('showcurrent_Callback',gcbo,[],guidata(gcbo)); end
        if ishandle(2), 
            dc2dinvres('compare_Callback',gcbo,[],guidata(gcbo)); 
        end
        if malstat==14, dc2dinvres('showip_Callback',gcbo,[],guidata(gcbo)); end
    end
end


% --------------------------------------------------------------------
function varargout = malfeld_Model(h, eventdata, handles, varargin)
global Mod N malstat FIX XX ZZ CHIQ RMS INV
cp=get(handles.malfeld,'CurrentPoint');
if isfield(N,'topo'),
    cp(1,2)=interp1(N.topo(:,1),N.topo(:,2),cp(1,1),'linear','extrap')-cp(1,2);
end
nx=max(find(Mod.x<cp(1,1)));
nz=max(find(Mod.z<cp(1,2)));
if (nx>=1)&&(nx<length(Mod.x))&&(nz>=1)&&(nz<length(Mod.z)), %found
    if malstat==11, % model cell resolution
        nr=(nz-1)*size(Mod.M,1)+nx;
        set(handles.slider,'Value',nr);
        dc2dinvres('slider_Callback',gcbo,[],guidata(gcbo));
    else
        st=get(handles.widerstand,'String');
        nu=str2double(st);
        stip=get(handles.ipmodelinput,'String');
        nuip=str2double(stip);
        if isnan(nu), nu=0; end
        if (nu>0)||(nuip>0)||ismember(st,'FUCXZ'),
            finalRect=rbbox;
            cp1=get(handles.malfeld,'CurrentPoint');
            if isfield(N,'topo'),
                cp1(1,2)=interp1(N.topo(:,1),N.topo(:,2),cp1(1,1),'linear','extrap')-cp1(1,2);
            end
            startx=max(find(Mod.x<cp(1,1)));
            startz=max(find(Mod.z<cp(1,2)));
            nx=max(find(Mod.x<cp1(1,1)));
            nz=max(find(Mod.z<cp1(1,2)));
            if startx>nx, du=startx;startx=nx;nx=du; end
            if startz>nz, du=startz;startz=nz;nz=du; end
            if nx>size(Mod.M,1), nx=size(Mod.M,1); end
            if nz>size(Mod.M,2), nz=size(Mod.M,2); end
            if nu>0,
                Mod.M(startx:nx,startz:nz)=nu;
                un=unique(Mod.M(:));
                if length(un)==1, % homogeneous case
                    Mod.R(:)=un;Mod.Lay(:)=un; 
                    RMS=rms(N.r,Mod.R,INV.lolo);
                    CHIQ=chi2(N.r,Mod.R,N.err,INV.lolo);
                    Mod.isfor=1;
                    dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo));
                else Mod.isfor=0; end
            else
                if ~isequal(size(Mod.M),size(FIX)), FIX=zeros(size(Mod.M)); end
                if ~isequal(size(Mod.M),size(XX)+[1 0]), XX=zeros(size(Mod.M)-[1 0]); end
                if ~isequal(size(Mod.M),size(ZZ)+[0 1]), ZZ=zeros(size(Mod.M)-[0 1]); end
                nr=0;
                if isequal(st,'F'), nr=-1; end
                if isequal(st,'C'), % compound model cells
                   nr=1; 
                   while(find(FIX==nr)), nr=nr+1; end
                end
                if isequal(st,'X'),
                    XX(startx:nx-1,startz:nz)=1;
                elseif isequal(st,'Z'),
                    ZZ(startx:nx,startz:nz-1)=1;
                else
                    FIX(startx:nx,startz:nz)=nr;
                end
            end
            if (nuip>0)&~ismember(st,'FUCXZ'), 
                if ~isfield(Mod,'IP')||(~isequal(size(Mod.M),size(Mod.IP))), 
                    Mod.IP=zeros(size(Mod.M)); end
                Mod.IP(startx:nx,startz:nz)=nuip; 
            end
            dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
%             set(handles.widerstand,'String',num2str(nu));
            save('model.mat','Mod');
            set(handles.backtoold,'Enable','On');    
        else
            if isfield(Mod,'IP')&&isequal(size(Mod.M),size(Mod.IP)),
                msgbox(sprintf('rho=%.1f ip=%.2f (Mod.x=%.1f-%.1f, Mod.z=%.1f-%.1f)\n',Mod.M(nx,nz),Mod.IP(nx,nz),Mod.x(nx),Mod.x(nx+1),Mod.z(nz),Mod.z(nz+1)),...
                    'Resistivity cell');
            else
                msgbox(sprintf('rho=%.1f (Mod.x=%.1f-%.1f, Mod.z=%.1f-%.1f)\n',Mod.M(nx,nz),Mod.x(nx),Mod.x(nx+1),Mod.z(nz),Mod.z(nz+1)),...
                    'Resistivity cell');
            end
        end    
    end
end

% --------------------------------------------------------------------
function varargout = slider_Callback(h, eventdata, handles, varargin)
global MAL malstat S Mod N RM RD INV
mal=MAL;
mal.log=0;mal.cmap=2;mal.cauto=0;mal.elec=0;
if malstat>0,
    nr=round(get(handles.slider,'Value'));
    set(handles.editslider,'String',num2str(nr));
    set(handles.slider,'Value',nr);
    axes(handles.malfeld);
end
nm=numel(Mod.M);
switch malstat,
    case 0, % model - lambda
        val=get(handles.slider,'Value');
        lam=10^(3*(1-val));
        if lam>100, ss=num2str(round(lam)); else ss=num2str(lam,'%.1f'); end
        set(handles.editslider,'String',ss);
%         dc2dinvres('editslider_Callback',gcbo,[],guidata(gcbo))
    case 9, % Sensitivity
        snr=Mod.M;
        if size(S,2)==nm, snr(:)=S(nr,:)'; end
        if size(S,1)==nm, snr(:)=S(:,nr); end
%         mm=max(abs(snr(:)))/10;
        mm=std(snr(:))*3;
        mal.cmin=-mm;mal.cmax=mm;
%         draw2dmodel(Mod.x,Mod.z,snr,mal);
        patch2dmodel(Mod.x,Mod.z,snr,mal,N);
        %     hold on
        %     plot(N.elec(N.a(nr),1),N.elec(N.a(nr),2),'kv');
        %     plot(N.elec(N.m(nr),1),N.elec(N.m(nr),2),'kv');    
        %     hold off
        ss=sprintf('A=%.1f',N.elec(N.a(nr),1));
        if N.b(nr)>0,
            ss=sprintf('%s B=%.1f',ss,N.elec(N.b(nr),1));
        end
        ss=sprintf('%s Mod.M=%.1f',ss,N.elec(N.m(nr),1));
        if N.n(nr)>0,
            ss=sprintf('%s N=%.1f',ss,N.elec(N.n(nr),1));
        end
        ss=''; %!!! (too big)
        ss=sprintf('%s rhoa=%.1f',ss,N.r(nr));
        if isfield(N,'err')&&(length(N.err)>=nr),
            ss=sprintf('%s Err=%.1f%%',ss,N.err(nr)*100);
        end
        set(handles.helptext,'String',ss);
    case 10, % Data Covariance
        [mids,seps,ii,kk]=showdata2d(N,RD(:,nr));
        hold on;plot(mids(ii(nr)),kk(nr),'kx','MarkerSize',10); hold off
    case 11, % Model Covariance
        if isempty(RM), 
            C=smooth2d1st(Mod.x,Mod.z,1);D=spdiags(1./log(1+N.err),0,length(N.err),length(N.err));
            P=1;%if INV.redu>1, P=pmatrix2d(Mod.x,Mod.z); end
            set(gcf,'Pointer','watch');
            rmn=modelcellres(S,INV.lam,nr,C,D,P);
            set(gcf,'Pointer','arrow');
        else rmn=RM(:,nr); end
        mm=max(abs(rmn));mal.cmin=-mm;mal.cmax=mm;
%         draw2dmodel(Mod.x,Mod.z,rmn,mal);
        patch2dmodel(Mod.x,Mod.z,rmn,mal,N);
        lx=length(Mod.x)-1;
        nx=mod(nr,lx);if nx==0, nx=lx; end
        nz=(nr-nx)/lx+1;
        hold on
        line([Mod.x(nx) Mod.x(nx+1) Mod.x(nx+1) Mod.x(nx) Mod.x(nx)],[Mod.z(nz) Mod.z(nz) Mod.z(nz+1) Mod.z(nz+1) Mod.z(nz)],...
            'LineWidth',2,'Color','Black');
        hold off
        im=get(handles.malfeld,'Children');
        for i=1:length(im),
            if strcmp(get(im(i),'type'),'surface'),
                set(im(i),'ButtonDownFcn','dc2dinvres(''malfeld_Model'',gcbo,[],guidata(gcbo))')
            end
        end
    otherwise
        set(handles.preview,'Visible','Off');
        set(handles.slider,'Visible','Off');
        set(handles.editslider,'Visible','Off');
        set(handles.dispslider,'Visible','Off');
end
if (malstat==9)|(malstat==10), % draw electrodes
    el=[N.elec(N.a(nr),:);N.elec(N.m(nr),:)];
    ch='AM';
    if N.b(nr)>0, el=[el;N.elec(N.b(nr),:)];ch=[ch 'B']; end
    if N.n(nr)>0, el=[el;N.elec(N.n(nr),:)];ch=[ch 'N']; end
    if malstat==10, el(:,2)=0.5; end
    xt=get(gca,'XTick');xl=cellstr(get(gca,'XTickLabel'));
    for k=1:length(ch),
        l=find(xt==el(k,1));
        if isempty(l), l=length(xt)+1; end
        xl{l}=ch(k);xt(l)=el(k,1);
    end
    [xt,l]=sort(xt);
    yl=xl;for i=1:length(xl), yl{i}=xl{l(i)}; end
    set(handles.malfeld,'XTick',xt,'XTickLabel',yl);
    hold on
    plot(el(:,1),el(:,2),'kv','MarkerSize',5);
    hold off
end % slider_callback    

% --------------------------------------------------------------------
function varargout = editslider_Callback(h, eventdata, handles, varargin)
global malstat INV
sl=get(handles.editslider,'String');
if malstat==0,
    lam=str2double(sl);INV.lam=lam; % set regularization parameter
    if lam<1, lam=1; end % only for the slider value!
    if lam>1000, lam=1000; end
    set(handles.slider,'Value',1-log10(lam)/3);
else
    nr=round(str2double(sl));
    if (nr>get(handles.slider,'Max'))||(nr<get(handles.slider,'Min')),
        nr=1;
    end
    set(handles.slider,'Value',nr);
    dc2dinvres('slider_Callback',gcbo,[],guidata(gcbo))
end

% --------------------------------------------------------------------
function varargout = widerstand_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function noisify_Callback(hObject, eventdata, handles)
global N Mod
daterr;
noise=randn(size(N.r)).*N.err;
Mod.R=Mod.R.*(1+noise);
dc2dinvres('showfor_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function seterror_Callback(hObject, eventdata, handles)
global N Mod CHIQ INV
daterr;
CHIQ=[CHIQ chi2(N.r,Mod.R,N.err,INV.lolo)];
messg(sprintf('Error estimate: min(err)=%.1f%% max(err)=%.1f%% chi^2=%.2f',...
    min(N.err)*100,max(N.err)*100,CHIQ(end)));

function loadsens_Callback(hObject,eventdata,handles)
global datfile S Mod
[fpath,name,ext]=fileparts(datfile);
matfile=strrep(datfile,ext,'-sens.mat');
if exist(matfile,'file'), 
    load(matfile);
    Mod.Cov=zeros(size(Mod.M));
    nm=numel(Mod.M);
    if size(S,2)==nm,
        %     Mod.Cov(:)=sum(abs(S));
        for i=1:nm, Mod.Cov(i)=sum(abs(S(:,i))); end
    end
    if size(S,1)==nm,
        %     Mod.Cov(:)=sum(abs(S),2);
        for i=1:nm, Mod.Cov(i)=sum(abs(S(i,:))); end
    end
    messg('Loading Sensitivity matrix from disk');
    if ~exist('xsave','var'), xsave=Mod.x;zsave=Mod.z; end % !!! Ausnahme
    if (~isequal(Mod.x(:),xsave(:)))||(~isequal(Mod.z(:),zsave(:))),
        S=[];Mod.Cov=[];
        messg('Recomputation of sensitivity will be necessary');
    end
end


% --------------------------------------------------------------------
function createdata_Callback(hObject, eventdata, handles)
global N Mod INV RMS datfile FOR CHIQ
pname=fileparts(datfile);
uiwait(ncreate);
datfile=fullfile(pname,'new.dat');
if isfield(N,'topo'), N=rmfield(N,'topo'); end
dc2dinvres('showdata_Callback',gcbo,[],guidata(gcbo));
if ~isfield(Mod,'M')||isempty(Mod.M)||...
        strcmp(questdlg('Create new model or keep old','Create new model','New','Keep','New'),'New')
    [Mod.x,Mod.z,Mod.M]=modelfromdata2d(N);
    rho=Mod.M(1,1);Mod.Lay=rho;
    Mod.R=ones(size(N.r))*rho;
else
    set(gcf,'Pointer','watch');
    Mod.R=dcfwd2d(Mod.x,Mod.z,Mod.M,Mod.Lay,N,FOR);
    set(gcf,'Pointer','arrow');
    dc2dinvres('showfor_Callback',gcbo,[],guidata(gcbo));
end
RMS=rms(N.r,Mod.R,INV.lolo);
CHIQ=chi2(N.r,Mod.R,N.err,INV.lolo);
dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
try 
    set(hObject,'Name',[ 'DC2DInvRes - ' datfile ]);
catch
    set(handles.inv2d,'Name',[ 'DC2DInvRes - ' datfile ]);
end


% --------------------------------------------------------------------
function comparemod_Callback(hObject, eventdata, handles)
global Mod datfile INV MAL FIX XX ZZ N
[fpath,name,ext]=fileparts(datfile);
if isempty(ext),
    outfile='*.mod';
else
    outfile=strrep(datfile,ext,'.mod');
end
[fname,pname]=uigetfile('*.mod','Model to compare',outfile);
if fname==0, return; end
[M1,x1,z1]=modelimport2d(fullfile(pname,fname));
if ~(isequal(Mod.x(:),x1(:))&&isequal(Mod.z(:),z1(:))),
    errordlg('Model parameterization is not identical!');
    return;
end
dM=M1./Mod.M;
axes(handles.malfeld);
mm=max(max(dM(:)),1/min(dM(:)));
mal=MAL;mal.cauto=0;mal.clog=1;mal.cmin=1./mm;mal.cmax=mm;
patch2dmodel(Mod.x,Mod.z,dM,mal);
hold on
if isequal(size(FIX),size(Mod.M)),
    [i,j]=find(FIX==-1);
    plot((Mod.x(i)+Mod.x(i+1))/2,(Mod.z(j)+Mod.z(j+1))/2,'wx');
    [i,j]=find(FIX>0);
    plot((Mod.x(i)+Mod.x(i+1))/2,(Mod.z(j)+Mod.z(j+1))/2,'wo');
end
if isequal(size(XX)+[1 0],size(Mod.M)),
    [i,j]=find(XX>0);
    plot(Mod.x(i+1),(Mod.z(j)+Mod.z(j+1))/2,'w+');
end
if isequal(size(ZZ)+[0 1],size(Mod.M)),
    [i,j]=find(ZZ>0);
    plot((Mod.x(i)+Mod.x(i+1))/2,Mod.z(j+1),'w+');    
end
hold off
if INV.lolo,
    ss=sprintf('Model Difference RMS = %.2f%%',rms(log(Mod.M(:)),log(M1(:))));
else
    ss=sprintf('Model Difference RMS = %.2f%%',rms(Mod.M(:),M1(:)));
end
% uiwait(msgbox(ss,'Compare models'));
if strcmp(questdlg('Export ratio as model file?',ss),'Yes'),
    [fpath,name,ext]=fileparts(datfile);
    [fname,pname]=uiputfile(fullfile(fpath,'ratio.mod'),'Export model to ASCII-File');
    if fname~=0,
        outfile=fullfile(pname,fname);
        [fpath,name,ee]=fileparts(outfile);
        if strcmp(ee,''), outfile=[outfile '.mod']; end
        sec=Mod.Cov;isip=0;
        if isfield(Mod,'IP')&&isequal(size(Mod.M),size(Mod.IP)), 
            sec=Mod.IP;isip=1; end
        if isfield(N,'topo'),
            xz=interp1(N.topo(:,1),N.topo(:,2),Mod.x,'linear','extrap');
        else
            xz=zeros(size(Mod.x));
        end
        modelexport2d(outfile,1./dM,Mod.x,Mod.z,sec,0,xz,isip);
    end    
end


% --------------------------------------------------------------------
function importmod_Callback(hObject, eventdata, handles)
global Mod datfile S N RMS CHIQ FOR INV
[fpath,name,ext]=fileparts(datfile);
if isempty(ext),
    outfile='*.mod';
else
    outfile=strrep(datfile,ext,'.mod');
end
[fname,pname]=uigetfile(outfile,'Import model from ASCII-File');
if fname==0, return; end
[Mod.M,Mod.x,Mod.z]=modelimport2d(fullfile(pname,fname));
messg(['Importing model file ' fullfile(pname,fname)]);
Mod.Lay=ones(length(Mod.z),1);
for k=1:size(Mod.M,2),
    Mod.Lay(k)=median(Mod.M(:,k));
end
Mod.Lay(k+1:end)=median(Mod.M(:,end));
%Mod.Mref=Mod.M; % keine gute Idee, oder?
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
oldR=Mod.R;
if isstruct(N),
    dc2dinvres('forwarding_Callback',gcbo,[],guidata(gcbo));
%     set(gcf,'Pointer','watch');drawnow;
%     Mod.R=dcfwd2d(Mod.x,Mod.z,Mod.M,Mod.Lay,N,FOR);
%     RMS=[RMS rms(N.r,Mod.R,INV.lolo)];
%     CHIQ=[CHIQ chi2(N.r,Mod.R,N.err,INV.lolo)];
%     dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo));
%     if ishandle(2), dc2dinvres('compare_Callback',gcbo,[],guidata(gcbo)); end
    rho=median(N.r);dM=log(Mod.M(:))-log(rho);
    br=1;
    if ~isequal(sort(size(S)),sort([length(N.r) numel(Mod.M)])),
        messg('Sensitivity recalculation necessary!');
        answer=questdlg('Sensitivity wrong! Compute sensitivity now?');
        if strcmp(answer,'Yes'),
            dc2dinvres('compsens_Callback',gcbo,[],guidata(gcbo));
        else
            br=0;
        end
    end
    if br,
        % S = S + ((log(Mod.R)-log(rho)-S*dM)*(dM'))/(dM'*dM); % broyden
        if (size(S,1)==length(Mod.R))&&(size(S,2)==length(dM)),
            ddr=(log(Mod.R)-log(rho)-S*dM(:))/(dM(:)'*dM(:));
            for i=1:size(S,2), S(:,i)=S(:,i)+ddr*dM(i); end
            messg('Doing Broyden update of sensitivity matrix!');
        end
        if (size(S,2)==length(Mod.R))&&(size(S,1)==length(dM)), %trans+sparse
            ddr=(log(Mod.R)-log(rho)-dM(:)'*S)/(dM(:)'*dM(:));
            for i=1:size(S,1), S(i,:)=S(i,:)+ddr*dM(i); end
            messg('Doing Broyden update of sensitivity matrix!');
        end
    end
end
set(gcf,'Pointer','arrow');drawnow;

function zeile=destrip(zeile)
% strip string from comments (with # character)
aa=strfind(zeile,'#');
if ~isempty(aa), zeile=zeile(1:aa(1)-1); end


% --------------------------------------------------------------------
function editdata_Callback(hObject, eventdata, handles)
global N Mod RMS INV CHIQ
le=length(N.r);
ndelete;
N=deldeadelecs(N);
if length(N.r)~=le,
    RMS=[RMS rms(N.r,Mod.R,INV.lolo)];
    CHIQ=[CHIQ chi2(N.r,Mod.R,N.err,INV.lolo)];
    dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo));
end
dc2dinvres('showdata_Callback',gcbo,[],guidata(gcbo));


% --------------------------------------------------------------------
function help_Callback(hObject, eventdata, handles)
ibrowse(['file://' pwd filesep 'doc' filesep 'index.html']);

% --------------------------------------------------------------------
function rmsout_Callback(hObject, eventdata, handles)
global RMS CHIQ
messg(['CHI^2 = ',sprintf('%.2f ',CHIQ) '(RMS=' sprintf('%.2f%%',RMS(end)) ')']);


% --------------------------------------------------------------------
function backtoold_Callback(hObject, eventdata, handles)
global CHIQ RMS% Mod
if exist('model.mat','file'),
    load model.mat
    if ~isempty(CHIQ), CHIQ(end)=[]; end
    if ~isempty(RMS), RMS(end)=[]; end
    dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo));
    dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
else
    messg('No old model available');
end
set(handles.backtoold,'Enable','Off');

% --------------------------------------------------------------------
function referencemod_Callback(hObject, eventdata, handles)
global Mod
Mod.Mref=Mod.M;

% --------------------------------------------------------------------
function onlinehelp_Callback(hObject, eventdata, handles)
ibrowse('http://www.resistivity.net/dc2dinvres/');

% --------------------------------------------------------------------
function webinversion_Callback(hObject, eventdata, handles)
ibrowse('http://webinv.resistivity.net');

function blockify_Callback(hObject, eventdata, handles)
global Mod N FOR
save('model.mat','Mod');
set(handles.backtoold,'Enable','On');
% Mod.M=blockify2dmodel(Mod.M);
% Mod.M=cluster2dmodel(Mod.M);
data=log10(Mod.M(:));
isip=isfield(Mod,'IP')&&isequal(size(Mod.M),size(Mod.IP));
if isip, 
    data(:,2)=Mod.IP(:);%log10(IP(:)); 
end
%%
kruemm=0;mal=struct('clog',0,'alfa',1);
for i=1:15,
    param.c=i;
    result{i}=fcmcluster(data,param);    
    y(i)=result{i}.cost(end);
    if i>1,
        cl1=Mod.M;mf1=Mod.M;[mf1(:),cl1(:)]=max(result{i}.data.f,[],2);
        alfa1=(mf1*param.c-1)/(param.c-1);
        patch2dmodel(Mod.x,Mod.z,cl1,mal,N,alfa1);
        pause(1);
    end
    if i>2,
        ys=(y(3:end)-y(1:end-2))/2;
        kruemm(i-1)=(y(i)+y(i-2)-2*y(i-1))/((1+(y(i)-y(i-2))^2/4)^1.5);
        if kruemm(i-1)<kruemm(i-2), k=i-1;break; end
    end
end
snum=inputdlg('Please specify number of clusters (0=cancel)?','Cluster number',1,{num2str(k)});
k=str2double(snum{1});
if k==0, return; end
cl1=Mod.M;mf1=Mod.M;[mf1(:),cl1(:)]=max(result{k}.data.f,[],2);
alfa1=(mf1*param.c-1)/(param.c-1);
Mod.M(:)=10.^result{k}.cluster.v(cl1,1);
if isip, 
%     Mod.IP(:)=10.^result{k}.cluster.v(cl1,2); 
    Mod.IP(:)=result{k}.cluster.v(cl1,2); 
end
%%
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
dc2dinvres('forwarding_Callback',gcbo,[],guidata(gcbo));
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
if strcmp(questdlg('perform cluster postiteration?'),'Yes'),
    un=unique(Mod.M(:));
    fak=1.1;
    A=zeros(length(N.r),length(un));
    for i=1:length(un),
        M1=Mod.M;
        M1(Mod.M==un(i))=un(i)*fak;
        R1=dcfwd2d(Mod.x,Mod.z,M1,Mod.Lay,N,FOR);
        A(:,i)=(log(R1)-log(Mod.R))/log(fak);
    end
    dR=log(N.r)-log(Mod.R);
    dm=(A'*A+0*eye(length(un)))\(A'*dR);
    newM=exp(dm).*un;
    for i=1:length(un),
        Mod.M(Mod.M==un(i))=newM(i);
    end
    dc2dinvres('forwarding_Callback',gcbo,[],guidata(gcbo));
    dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
end

% --------------------------------------------------------------------
function showvoltage_Callback(hObject, eventdata, handles)
global N malstat
if isfield(N,'u')&&(length(N.u)==length(N.a)),
    showdata2d(N,abs(N.u)*1000,struct('log',1));
    set(handles.helptext,'Visible','On',...
        'String','Voltage in mV');
    malstat=12;
    im=get(handles.malfeld,'Children');
    set(im,'ButtonDownFcn','dc2dinvres(''malfeld_ButtonDownFcn'',gcbo,[],guidata(gcbo))')
    set(handles.preview,'Visible','Off');
    set(handles.slider,'Visible','Off');
    set(handles.editslider,'Visible','Off');
    set(handles.dispslider,'Visible','Off');
    set(handles.helptext,'String','U in mV - Click to see/delete data');
    set(handles.helptext,'Visible','On');
    set(handles.widerstand,'Visible','Off');
    set(handles.ipmodelinput,'Visible','Off');
end

% --------------------------------------------------------------------
function showcurrent_Callback(hObject, eventdata, handles)
global N malstat
if isfield(N,'i')&&(length(N.i)==length(N.a)),
    showdata2d(N,abs(N.i)*1000,struct('log',1));
    set(handles.helptext,'Visible','On',...
        'String','Current in mA');
    malstat=13;
    im=get(handles.malfeld,'Children');
    set(im,'ButtonDownFcn','dc2dinvres(''malfeld_ButtonDownFcn'',gcbo,[],guidata(gcbo))')
    set(handles.preview,'Visible','Off');
    set(handles.slider,'Visible','Off');
    set(handles.editslider,'Visible','Off');
    set(handles.dispslider,'Visible','Off');
    set(handles.helptext,'String','I in mA - Click to see/delete data');
    set(handles.helptext,'Visible','On');
    set(handles.widerstand,'Visible','Off');
    set(handles.ipmodelinput,'Visible','Off');
end


% --------------------------------------------------------------------
function ipinv_Callback(hObject, eventdata, handles)
global Mod INV S N
Mod.IP=Mod.M;
if isequal(sort(size(S)),sort([length(N.ip) numel(Mod.IP)])),
    if isfield(INV,'blocky')&&(INV.blocky>0)&&(INV.weight==1),
        dM0=log(Mod.M(:));
        if isfield(Mod,'Mref'), dM0=dM0-log(Mod.Mref(:)); end
        D=spdiags(1./log(1+N.err),0,length(N.err),length(N.err));
        rbzfak=1;if isfield(INV,'rbzfak'), rbzfak=INV.rbzfak; end
        rbnorm=1;if isfield(INV,'rbnorm'), rbnorm=INV.rbnorm; end
        [C,Cx,Cz]=smooth2d1st(Mod.x,Mod.z,rbnorm,rbzfak);
        Cxz=[Cx;Cz];
        sm=abs(Cxz*dM0);su2=sum(sm.^2);sua=sum(sm);
        wxz=su2/sua./sm;
        Cxz=spdiags(wxz,0,length(wxz),length(wxz))*Cxz;
        C=Cxz'*Cxz;
        if length(N.ip)==size(S,2),
            Mod.IP(:)=cglscdpt(S,N.ip,INV.lam,C,D);
        else
            Mod.IP(:)=cglscdp(S,N.ip,INV.lam,C,D);
        end
    else
        ipinv=INV;%inv.auto=2;
        fi=find(N.ip>0);
        if length(N.ip)==size(S,2), % transposed
            Mod.IP(:)=invers(S(:,fi),N.ip(fi),ipinv);
        else
%             Mod.IP(:)=exp(invers(S(fi,:),log(N.ip(fi)),ipinv));
            Mod.IP(:)=invers(S(fi,:),N.ip(fi),ipinv);
            % Mod.IP(:)=invers(S,N.ip,ipinv);
        end
        Mod.IP(Mod.IP<0)=0;
    end
    if size(S,2)==numel(Mod.IP), Mod.ipfor=S*Mod.IP(:); end
    if size(S,1)==numel(Mod.IP), Mod.ipfor=(Mod.IP(:)'*S)'; end
    ipchi=chi2(N.ip+eps,Mod.ipfor,N.err,0);
%     iprms=rms(N.ip+eps,Mod.ipfor,0);
    iprms=sqrt(mean((N.ip-Mod.ipfor).^2));
    messg(sprintf('IP misfit: chi^2=%.1f rms=%.1fmrad',ipchi,iprms));
    dc2dinvres('showcharge_Callback',gcbo,[],guidata(gcbo))
else
    errordlg('Sensitivity matrix size mismatching');
end

% --------------------------------------------------------------------
function showip_Callback(hObject, eventdata, handles)
global N malstat MAL
axes(handles.malfeld);
if isfield(N,'ip'),
    mal=MAL;mal.cauto=1;mal.canot='IP/mrad';mal.clog=0;mal.log=0;
    showdata2d(N,N.ip,mal);
%     set(handles.helptext,'Visible','On',...
%         'String','IP-Data in mrad');
    malstat=14;
    im=get(handles.malfeld,'Children');
    set(im,'ButtonDownFcn','dc2dinvres(''malfeld_ButtonDownFcn'',gcbo,[],guidata(gcbo))')
    set(handles.preview,'Visible','Off');
    set(handles.slider,'Visible','Off');
    set(handles.editslider,'Visible','Off');
    set(handles.dispslider,'Visible','Off');
    set(handles.helptext,'String','IP in mrad - Click to see/delete data');
    set(handles.helptext,'Visible','On');
    set(handles.widerstand,'Visible','Off');
    set(handles.ipmodelinput,'Visible','Off');
end

% --------------------------------------------------------------------
function showcharge_Callback(hObject, eventdata, handles)
global Mod N MAL malstat FIX XX ZZ
axes(handles.malfeld);
if isfield(Mod,'IP')&&isequal(size(Mod.IP),size(Mod.M)),
    mal=MAL;mal.cauto=1;mal.clog=(min(Mod.IP(:))>0);
    mal.log=mal.clog;mal.canot='phi/mrad';
    if isfield(MAL,'alpha')&&(MAL.alpha)&&isfield(Mod,'Cov'), 
        patch2dmodel(Mod.x,Mod.z,Mod.IP,mal,N,Mod.Cov);
    else patch2dmodel(Mod.x,Mod.z,Mod.IP,mal,N); end
    plotconstraints(Mod,FIX,XX,ZZ,N);    
    malstat=15;
    im=get(handles.malfeld,'Children');
    set(im,'ButtonDownFcn','dc2dinvres(''malfeld_Model'',gcbo,[],guidata(gcbo))')
    set(handles.preview,'Visible','Off');
    set(handles.slider,'Visible','Off');
    set(handles.editslider,'Visible','Off');
    set(handles.dispslider,'Visible','Off');
    set(handles.helptext,'String','IP model');
    set(handles.helptext,'Visible','On');
    set(handles.widerstand,'Visible','Off');
    set(handles.ipmodelinput,'Visible','Off');
end

% --------------------------------------------------------------------
function data_Callback(hObject, eventdata, handles)
global N Mod
sdata=sonoff(isstruct(N));
set(handles.seterror,'Enable',sdata);
set(handles.reweight,'Enable',sdata);
set(handles.editdata,'Enable',sdata);
set(handles.filewrite,'Enable',sdata);
set(handles.savelb,'Enable',sdata);
set(handles.exportohm,'Enable',sdata);
sdata=sonoff(isfield(N,'a')&&(length(N.a)==length(Mod.R)));
set(handles.setasdata,'Enable',sdata);
set(handles.noisify,'Enable',sdata);
sdata=sonoff(isfield(N,'a')&&isfield(Mod,'M')&&~isempty(Mod.M));
set(handles.forwarding,'Enable',sdata);

% --------------------------------------------------------------------
function compip_Callback(hObject, eventdata, handles)
global N MAL S Mod
if ~isfield(Mod,'ipfor'),
    Mod.ipfor=zeros(size(N.a));
    if size(S,2)==numel(Mod.IP), Mod.ipfor=S*Mod.IP(:); end
    if size(S,1)==numel(Mod.IP), Mod.ipfor=(Mod.IP(:)'*S)'; end
end
dip=N.ip-Mod.ipfor;
% dip=(1-Mod.ipfor(:)./N.ip(:))*100;
figure(2);
set(2,'Menubar','none','NumberTitle','off','Name','IP Data,Model response & misfit');
iconify(2);
mal=MAL;mal.cauto=0;
mal.cmax=max([N.ip(:);Mod.ipfor(:)]);
mal.cmin=min([N.ip(:);Mod.ipfor(:)]);
if mal.cmin<=0, 
    mal.clog=0;
    mal.cmax=max(abs([mal.cmin mal.cmax]));
    mal.cmin=-mal.cmax; 
end
subplot(3,1,1);
showdata2d(N,N.ip,mal);
xl=get(gca,'XLim');xl=xl(1)+0.05*diff(xl);
yl=get(gca,'YLim');yl=yl(1)+0.9*diff(yl);
text(xl,yl,'Measured IP Data');
subplot(3,1,2);
showdata2d(N,Mod.ipfor,mal);
text(xl,yl,'Calculated IP Data');
subplot(3,1,3);
mal.log=0;mal.cauto=0;
mm=max(abs(dip));mal.cmin=-mm;mal.cmax=mm;
showdata2d(N,dip,mal);
%colormap(jet);
text(xl,yl,'Misfit in mrad');


% --------------------------------------------------------------------
function modelfix_Callback(hObject, eventdata, handles)
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
set(handles.widerstand,'String','F');

% --------------------------------------------------------------------
function modelunfix_Callback(hObject, eventdata, handles)
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
set(handles.widerstand,'String','U');


% --------------------------------------------------------------------
function reweight_Callback(hObject, eventdata, handles)
global N Mod CHIQ INV RMS
dR=(log(N.r)-log(Mod.R))./log(1+N.err);
w=abs(dR)*sum(abs(dR))/sum(dR.^2);
w(w<1)=1;
N.err=N.err.*w;
RMS=[RMS rms(N.r,Mod.R,INV.lolo)];
CHIQ=[CHIQ chi2(N.r,Mod.R,N.err,INV.lolo)];
dc2dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
showdata2d(N,w);


% --- Executes during object creation, after setting all properties.
function ipmodelinput_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function ipmodelinput_Callback(hObject, eventdata, handles)


% --- Executes on button press in preview.
function preview_Callback(hObject, eventdata, handles)
global MAL INV Mod S N %Mod.Cov
INV.lam=str2double(get(handles.editslider,'String'));
dR=log(N.r)-log(Mod.R);
dM=invers(S,dR,INV,N.err);
nM=Mod.M;nM(:)=nM(:).*exp(dM);
axes(handles.malfeld);
patch2dmodel(Mod.x,Mod.z,nM,MAL,N);

function previewon(handles)
global INV
return
set(handles.preview,'Visible','On');
set(handles.slider,'Visible','On','min',0,'max',1,'Value',1-log10(INV.lam)/3);
set(handles.editslider,'Visible','On','String','1');
set(handles.dispslider,'Visible','On','String','smooth model ---- structured model');
dc2dinvres('slider_Callback',gcbo,[],guidata(gcbo))

function previewoff(handles)
set(handles.preview,'Visible','Off');
set(handles.editslider,'Visible','Off');
set(handles.dispslider,'Visible','Off');
dc2dinvres('slider_Callback',gcbo,[],guidata(gcbo))


% --------------------------------------------------------------------
function exportohm_Callback(hObject, eventdata, handles)
global datfile N
[pname,fname,ext]=fileparts(datfile);
outfile=strrep(datfile,ext,'.ohm');
[fname,pname]=uiputfile(outfile,'Save Ohm-File');
if fname~=0,
    ohmfile=strrep(fullfile(pname,fname),[pwd filesep],'');
    [pname,fname,ee]=fileparts(ohmfile);
    if strcmp(ee,''), ohmfile=[ohmfile '.ohm']; end
    messg(['Exporting to ohm-xz file ',ohmfile]);
    if isfield(N,'topo'),
        T=rmfield(N,'topo');
        [elec,ii]=sort(N.elec(:,1));
        jj=1:length(ii);jj(ii)=1:length(ii);
        fa=find(N.a);T.a(fa)=jj(N.a(fa));
        fb=find(N.b);T.b(fb)=jj(N.b(fb));
        fm=find(N.m);T.m(fm)=jj(N.m(fm));
        fn=find(N.n);T.n(fn)=jj(N.n(fm));
        typ=questdlg('Is the topography measured on tape or real coordinates(xz)?','Topography type','Tape','XZ','Tape');
        if isequal(typ,'Tape'),
            %         T.elec=N.elec;
            %         T.elec(:,2)=interp1(N.topo(:,1),N.topo(:,2),N.elec(:,1),'linear','extrap');;
            T.elec=mbm2xz(elec(:,1),N.topo,0);
        else % replaced N.elec by (sorted) elec
            T.elec=mbm2xz(elec(:,1),N.topo,1);
        end
    else
        T=N;
    end
    T.rho=N.r./N.k;
    if isfield(T,'r'), T=rmfield(T,'r'); end
    if isfield(T,'i'), T=rmfield(T,'i'); end
    if isfield(T,'u'), T=rmfield(T,'u'); end
%     if isfield(T,'ip'), T=rmfield(T,'ip'); end
    saveinv2dfile(ohmfile,T);    
end

% --------------------------------------------------------------------
function vtkexport_Callback(hObject, eventdata, handles)
global N Mod datfile
[path,name,ext]=fileparts(datfile);
outfile=strrep(datfile,ext,'.vtk');
[fname,pname]=uiputfile(outfile,'Export model to VTK-File');
if fname==0, return; end% no valid name
outfile=fullfile(pname,fname);
[ff,pp,ee]=fileparts(outfile);
if strcmp(ee,''), outfile=[outfile '.vtk']; end
messg(['Exporting model to VTK File : ' outfile]);
if isfield(N,'topo')&&(size(N.topo,2)>1),
    vtkexport2d(outfile,Mod.M,Mod.x,Mod.z,N.topo);
else
    vtkexport2d(outfile,Mod.M,Mod.x,Mod.z);
end


% --------------------------------------------------------------------
function timelapse_Callback(hObject, eventdata, handles)
global Mod INV datfile MAL S N
[pname,fname,ext]=fileparts(datfile);
infile=strrep(datfile,[fname ext],'*.dat');
[newfile,newpath]=uigetfile(infile,'Time lapse file!');
N1=read2dfile(fullfile(newpath,newfile));
dM=Mod.M;
if ~isequal(N.elec,N1.elec), uiwait(errordlg('Electrodes not identical!','Time lapse error'));return; end
[c,ia,ib]=intersect([N1.a N1.b N1.m N1.n],[N.a N.b N.m N.n],'rows');
if ~any(ia), return; end
if INV.lolo, dR=log(N1.r(ia))-log(N.r(ib));
else dR=N1.r(ia)-N.r(ib); end
if size(S,1)==length(N.r), dM(:)=exp(invers(S(ib,:),dR,INV,(N1.err(ia)+N.err(ib))/2)); end
if size(S,2)==length(N.r), dM(:)=exp(invers(S(:,ib),dR,INV,(N1.err(ia)+N.err(ib))/2)); end
ggN.elec=N.elec;ggN.a=N.a(ib);ggN.b=N.b(ib);ggN.m=N.m(ib);
ggN.n=N.n(ib);ggN.err=(N1.err(ia)+N.err(ib))/2;ggN.r=N1.r(ia)./N.r(ib);
saveinv2dfile(fullfile(newpath,'ratio.dat'),ggN);
% only for identical configurations
% if length(N1.r)~=length(N.r), return; end
% if INV.lolo, dR=log(N1.r(:))-log(N.r(:));
% else dR=N1.r(:)-N.r(:); end
% dM(:)=exp(invers(S,dR,INV,N.err));
mal=MAL;mal.cauto=0;
if INV.lolo,
    mm=max(max(dM(:)),1/min(dM(:)));
    mal.clog=1;mal.cmin=1./mm;mal.cmax=mm;
else
    mm=max(abs(dM));
    mal.clog=0;mal.cmin=-mm;mal.cmax=mm;
end
patch2dmodel(Mod.x,Mod.z,dM,mal);
vtkexport2d('timelapse.vtk',dM,Mod.x,Mod.z);

% --------------------------------------------------------------------
function constraints_Callback(hObject, eventdata, handles)
set(handles.widerstand,'String','0');
global Mod
smod=sonoff((~isfield(Mod,'M'))||(~isempty(Mod.M)));
set(handles.decouplex,'Enable',smod);
set(handles.decouplez,'Enable',smod);
set(handles.referencemod,'Enable',smod);
set(handles.borefile,'Enable',smod);
set(handles.readline,'Enable',smod);
set(handles.removedec,'Enable',smod);
set(handles.removefixes,'Enable',smod);
set(handles.modelfix,'Enable',smod);
set(handles.modelunfix,'Enable',smod);
set(handles.compoundcells,'Enable',smod);


% --------------------------------------------------------------------
function decouplex_Callback(hObject, eventdata, handles)
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
set(handles.widerstand,'String','X');


% --------------------------------------------------------------------
function decouplez_Callback(hObject, eventdata, handles)
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
set(handles.widerstand,'String','Z');

% --------------------------------------------------------------------
function removedec_Callback(hObject, eventdata, handles)
global XX ZZ
XX(:)=0;ZZ(:)=0;
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function removefixes_Callback(hObject, eventdata, handles)
global FIX
FIX(:)=0;
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function borefile_Callback(hObject, eventdata, handles)
global datfile ZZ Mod
[pp,ff,ee]=fileparts(datfile);
infile=strrep(datfile,ee,'.bor');
[fname,pname]=uigetfile({'*.bor','BOR-File (*.bor)';'*.*','All Files (*.*)'},...
    'Choose borehole file',infile);
if ~ischar(fname), return; end
borfile=fullfile(pname,fname);
Bor=loadborfile(borfile);
for i=1:length(Bor.lay),
    Bor.nlay{i}=zeros(size(Bor.lay{i}));
    for j=1:length(Bor.lay{i}),
        [mi,k]=min(abs(Mod.z-Bor.lay{i}(j)));
        Bor.nlay{i}(j)=k(1);
    end
    Bor.nlay{i}(diff(Bor.nlay{i})==0)=[];
end
sm=size(Mod.M)-[0 1];
if ~isequal(size(ZZ),sm), ZZ=zeros(sm); end
for i=1:length(Bor.nlay),
    imi=max(find(Mod.x<Bor.pos(i)));
    ima=min(find(Mod.x>Bor.pos(i)));
    while(ima-imi<1), imi=imi-1;ima=ima+1; end
    for j=1:length(Bor.nlay{i}),
        while(imi>1)&&(ima<=sm(1))&&(Mod.x(ima)-Mod.x(imi)<Mod.z(Bor.nlay{i}(j))), 
            imi=imi-1;ima=ima+1; end
        if Bor.nlay{i}(j)>1, ZZ(imi:ima-1,Bor.nlay{i}(j)-1)=1; end
    end
end
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo))

function erg=sonoff(value)
erg='Off';
if (nargin>0)&&isnumeric(value)&&(value>0), erg='On'; end
if (nargin>0)&&islogical(value)&&(value), erg='On'; end


% --------------------------------------------------------------------
function compoundcells_Callback(hObject, eventdata, handles)
set(handles.widerstand,'String','C');
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));


% --------------------------------------------------------------------
function pickcell_Callback(hObject, eventdata, handles)
set(handles.widerstand,'String','0');
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));


% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
f=msgbox({'DC2dInvRes';'Version 2.12';'Author: Thomas Günther';'Resistivity.net'},'About');
iconify(f);uiwait(f);


% --------------------------------------------------------------------
function underwater_Callback(hObject, eventdata, handles)
global Mod N FIX ZZ
nmz=size(Mod.M,2);
if ~isequal(size(Mod.M),size(FIX)), FIX=zeros(size(Mod.M)); end
if ~isequal(size(Mod.M)-[0 1],size(ZZ)), ZZ=zeros(size(Mod.M)-[0 1]); end
xm=Mod.x(1:end-1)+diff(Mod.x)/2;
zz=interp1(N.elec(:,1),N.elec(:,2),xm,'nearest','extrap');
for i=1:length(zz),
   k=max(find(Mod.z<=zz(i)))-1;
   if k>0, 
       FIX(i,1:min(k,nmz))=-1; 
%        ZZ(i,min(k,nmz-1))=1;
%        ZZ(i,min(k,nmz-2)+(0:1))=1;
       dz=Mod.z(k+1)-Mod.z(k);
       if (zz(i)-Mod.z(k+1))/dz>0.3, ZZ(i,min(k+1,nmz-1))=1; end
       if (Mod.z(k+2)-zz(i))/dz>0.3, ZZ(i,min(k,nmz-1))=1; end
   end
end
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));


% --------------------------------------------------------------------
function rollalong_Callback(hObject, eventdata, handles)
global N Mod INV MAL FIX ZZ FOR
if isequal(size(Mod.M),size(FIX)), FFIX=FIX; else FFIX=[]; end
if isequal(size(Mod.M),size(ZZ)+[0 1]), ZZZ=ZZ; else ZZZ=[]; end
set(gcf,'Pointer','watch');
MAL.alpha=0;
xorg=Mod.x;FIX=[];ZZ=[];docov=0;FOR.silent=1;INV.silent=1;
ff=find(FFIX);fm=Mod.M(ff);
if ~isequal(size(Mod.Cov),size(Mod.M)), docov=1;Mod.Cov=zeros(size(Mod.M)); end
try,
    if min(diff(N.elec(:,1)))<0, 
        N=sort2delecs(N,1);% error('sort electrodes before'); 
    end
    if ~isfield(INV,'lbound'), INV.lbound=0; end
    if ~isfield(INV,'ubound'), INV.ubound=0; end
    
    abmn=[N.a N.b N.m N.n];
    ma=max(abmn,[],2);
    mi=min(abmn,[],2);
    dd=fix((max(ma-mi)+1))*2;
    del=median(diff(N.elec(:,1)));
    
    nel=size(N.elec,1);
%     rho=median(N.r);Mod.M(:)=rho;Mod.R=ones(size(N.r))*rho;
    patch2dmodel(xorg,Mod.z,Mod.M,MAL,N);
    l=0;mael=0;
    while mael<nel,
        miel=fix(l*dd/2+1);
        mael=min(miel+dd-1,nel);if nel-mael<15, mael=nel; end
        NN.elec=N.elec(miel:mael,:);
        fi=find((mi>=miel)&(ma<=mael));
        NN.a=N.a(fi)-miel+1;NN.b=N.b(fi)-miel+1;NN.m=N.m(fi)-miel+1;NN.n=N.n(fi)-miel+1;
        NN.r=N.r(fi);NN.k=N.k(fi);NN.err=N.err(fi);
        ixa=max(max(find(xorg<NN.elec(1,1)))-2,1);
        ixe=min(find(xorg>NN.elec(end,1)))+2;
        if isempty(ixe)||(ixe>length(xorg)), ixe=length(xorg); end
%         fprintf('\nElectrodes %d-%d (%.1f..%.1fm) ',miel,mael,xorg(ixa),xorg(ixe));
        messg(sprintf('Electrodes %d-%d (%.1f..%.1fm) ',miel,mael,xorg(ixa),xorg(ixe)));
        MM=Mod.M(ixa:ixe-1,:);Mod.x=xorg(ixa:ixe);nM=MM;
        MMref=Mod.Mref(ixa:ixe-1,:);
        %         RR=Mod.R(fi);   
        RR=abs(dcfwd2d(Mod.x,Mod.z,MM,Mod.Lay,NN,FOR));        
        if isempty(FFIX), FIX=zeros(size(MM)); else FIX=FFIX(ixa:ixe-1,:); end
        if isempty(ZZZ), ZZ=[]; else ZZ=ZZZ(ixa:ixe-1,:); end
        if l>0,
            fixto=fix(size(MM,1)/6);appendmessg(sprintf('fixed to %.1fm ',Mod.x(fixto+1)));
            FIX(1:fixto,:)=-1;
        end
%         fprintf('chi2=%.1f ',chi2(NN.r,RR,NN.err,INV.lolo));
        appendmessg(sprintf('chi2=%.1f',chi2(NN.r,RR,NN.err,INV.lolo)));drawnow;
        S=calcsens2d(Mod.x,Mod.z,NN,struct('silent',1));
%         Mod.Cov=zeros(size(Mod.M));Mod.Cov(:)=sum(abs(S));
        chiq=1e10;
        for i=1:6,
            dR=logtrans(NN.r,INV.lbound,INV.ubound)-logtrans(RR,INV.lbound,INV.ubound);
            dM0=logtrans(MM(:),INV.lbound,INV.ubound)-logtrans(MMref(:),INV.lbound,INV.ubound);
            dM=invers(S,dR,INV,NN.err,dM0);
            nM(:)=ilogtrans(logtrans(MM(:),INV.lbound,INV.ubound)+dM);
            R1=abs(dcfwd2d(Mod.x,Mod.z,nM,Mod.Lay,NN,FOR));
            fak=linesearch(NN,RR,R1,INV.lolo);
            if fak==0, 
                taug=0.3; % sampling point
                nM(:)=ilogtrans(logtrans(MM(:),INV.lbound,INV.ubound)+dM*taug,INV.lbound,INV.ubound);
                Ri=abs(dcfwd2d(Mod.x,Mod.z,nM,Mod.Lay,NN,FOR));
                chiqs=[chiq chi2(NN.r,Ri,NN.err,1) chi2(NN.r,R1,NN.err,1)];
                G=[1 0 0;1 taug taug^2;1 1 1];
                xx=G\chiqs';
                fak=-xx(2)/2/xx(3);
            end %appendmessg('->Line search failed.');break; end
            if fak<0.95,
                MM(:)=ilogtrans(logtrans(MM(:),INV.lbound,INV.ubound)+dM*fak);
                RR=abs(dcfwd2d(Mod.x,Mod.z,MM,Mod.Lay,NN,FOR));
            else
                MM=nM;RR=R1;
            end
            oldchiq=chiq;
            chiq=chi2(NN.r,RR,NN.err,INV.lolo);
            appendmessg(sprintf('->%.1f',chiq));pause(0.1);
            if (chiq<1)||(chiq/oldchiq>0.9), break; end
            if isfield(INV,'robust')&&(INV.robust>0),
                dR=(log(NN.r)-log(RR))./log(1+NN.err);
                w=abs(dR)*sum(abs(dR))/sum(dR.^2);w(w<1)=1;
                NN.err=NN.err.*w;
            end
        end
        appendmessg(sprintf('(%.1f%%)',rms(NN.r,RR,1)));
        Mod.M(ixa:ixe-1,:)=MM;Mod.R(fi)=RR;
        %% conductivity prolongation
        for k=1:size(MM,2), Mod.M(ixe:end,k)=MM(end,k); end
        Mod.M(ff)=fm;
        if docov, aa=Mod.Cov(ixa:ixe-1,:);aa(:)=sum(abs(S));
            Mod.Cov(ixa:ixe-1,:)=aa; end
        mal=MAL;mal.cauto=0;mal.cmin=interperc(MM(:),5);mal.cmax=interperc(MM(:),95);
%         plotlongmod(setfield(Mod,'x',xorg),mal);
        patch2dmodel(xorg,Mod.z,Mod.M,MAL,N);
        line(min(Mod.x)*[1 1],minmax(Mod.z),'Color','white');
        line(max(Mod.x)*[1 1],minmax(Mod.z),'Color','white')
        pause(0.1);        
        l=l+1;
    end
    Mod.x=xorg;FIX=FFIX;ZZ=ZZZ;FOR.silent=0;INV.silent=0;
    messg('successfully completed roll-along inversion');
catch,
    disp(lasterr);
    Mod.x=xorg;FIX=FFIX;ZZ=ZZZ;FOR.silent=0;INV.silent=0;
end
set(gcf,'Pointer','arrow');
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
set(figure(8),'MenuBar','none','NumberTitle','off','Name','LongProfile');
plotlongmod(Mod,MAL);


% --------------------------------------------------------------------
function readline_Callback(hObject, eventdata, handles)
global Mod ZZ datfile
[fpath,name,ext]=fileparts(datfile);
outfile=fullfile(fpath,'*.txt;*.csv');
[fname,pname]=uigetfile(outfile,'Read line from ascii file');
if ~ischar(fname), return; end
% A=load(fullfile(pname,fname));
% A(:,2)=abs(A(:,2));
fid=fopen(fullfile(pname,fname));
A=mytextscan(fid,'%f%f');
fclose(fid);
A{2}=abs(A{2});
nmz=size(Mod.M,2);
if ~isequal(size(Mod.M)-[0 1],size(ZZ)), ZZ=zeros(size(Mod.M)-[0 1]); end
xm=Mod.x(1:end-1)+diff(Mod.x)/2;
% zz=interp1(A(:,1),A(:,2),xm,'nearest','extrap');
zz=interp1(A{1},A{2},xm,'nearest');
for i=1:length(zz),
    if ~isnan(zz(i)),
        k=max(find(Mod.z<=zz(i)))-1;
        if k>0, 
            dz=Mod.z(k+1)-Mod.z(k);
            if (zz(i)-Mod.z(k+1))/dz>0.4, ZZ(i,min(k+1,nmz-1))=1; end
            if (Mod.z(k+2)-zz(i))/dz>0.4, ZZ(i,min(k,nmz-1))=1; end
        end
    end
end
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));


% --------------------------------------------------------------------
function releasenotes_Callback(hObject, eventdata, handles)
ibrowse('%23release.txt');

% --------------------------------------------------------------------
function readuline_Callback(hObject, eventdata, handles)
global Mod N FIX ZZ datfile
[fpath,name,ext]=fileparts(datfile);
outfile=fullfile(fpath,'*.txt;*.csv;*.xd');
[fname,pname]=uigetfile(outfile,'Read line from ascii file');
if ~ischar(fname), return; end
fid=fopen(fullfile(pname,fname));A=mytextscan(fid,'%f%f');fclose(fid);
A{2}=abs(A{2});
nmz=size(Mod.M,2);
if ~isequal(size(Mod.M),size(FIX)), FIX=zeros(size(Mod.M)); end
if ~isequal(size(Mod.M)-[0 1],size(ZZ)), ZZ=zeros(size(Mod.M)-[0 1]); end
xm=Mod.x(1:end-1)+diff(Mod.x)/2;
zz=interp1(A{1},A{2},xm,'nearest');
for i=1:length(zz),
    if ~isnan(zz(i)),
        k=max(find(Mod.z<=zz(i)))-1;
        if k>0, 
            dz=Mod.z(k+1)-Mod.z(k);
            if (zz(i)-Mod.z(k+1))/dz>0.4, ZZ(i,min(k+1,nmz-1))=1; end
            if (Mod.z(k+2)-zz(i))/dz>0.4, ZZ(i,min(k,nmz-1))=1; end
            FIX(i,1:min(k,nmz))=-1; 
        end
    end
end
dc2dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
