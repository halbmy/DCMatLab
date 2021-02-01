function varargout = dc3dinvres(varargin)
% dc3dinvres Application M-file for dc3dinvres.fig
%    FIG = dc3dinvres launch dc3dinvres GUI.
%    dc3dinvres('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 08-Jan-2007 13:18:57

if nargin == 0  % LAUNCH GUI
    
    fig = openfig(mfilename,'reuse');
  % Use system color scheme for figure:
    set(fig,'Color',get(0,'DefaultUicontrolBackgroundColor'));

  % Generate a structure of handles to pass to callbacks, and store it. 
    handles = guihandles(fig);
    guidata(fig, handles);
    iconify(fig);
    global Model S N MAL INV FOR output RMS datfile
    output=handles.output;
    verstr='2.11.6';
    set(handles.figure1,'Name',['DC3dInvRes v.' verstr ' - Thomas Günther']);
    set(handles.bg,'String','DC3dInvRes');
    mess=get(handles.output,'String');
    mess{2}=['Version ' verstr];
    set(handles.output,'String',mess);
    
    INV=struct('method',0,'weight',1,'redu',0,'mitschicht',0,'lolo',1,'lam',30,...
        'mico',0.4,'sens',0,'auto',2,'update',0,'const',1,'robust',0,'blocky',0,...
        'rbzfak',1,'glob',1,'linesearch',1,'lbound',0,'ubound',0,'spsens',0);
    FOR=struct('method',0,'acc',1e-4,'tol',1e-4,'maxit',50,'rand',4,...
        'prolong',4,'zusatz',2,'refine',1,'direct',-1,'fillup',1);
    MAL=struct('xy',1,'log',1,'cauto',1,'cmap',0,'cmin',25,'cmax',400,...
        'nu',0,'nv',0,'startwith',1,'elec',1,'xdir',0,'ydir',0,'cont',[],...
        'vie',0,'clog',1);
    if exist('setup.mat','file'), load('setup.mat'); end
    global libmmfile
    libmmfile=checklic3;
    if ~isequal(libmmfile,4),
        set(handles.save,'Enable','off');
        set(handles.modexp,'Enable','off');
        set(handles.exportvtk,'Enable','off');
        uiwait(warndlg({'No valid license file found! Some functions (export) are not available!',...
                'For a valid license file call Generate License Code from the start menu and send it to thomas@resistivity.net'}...
            ,'No valid license file')); 
    end 
    if nargout > 0
        varargout{1} = fig;
    end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
    try
        if (nargout)
            [varargout{1:nargout}] = feval(varargin{:});
        else
            feval(varargin{:});
        end
    catch
        disp(lasterr);
    end
end

% ---- INVERSION -----------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function reset_Callback(h,eventdata,handles,varargin)
global Model N MAL S RMS INV datfile CHIQ Coverage Mref
% Sigma_q=Mittelwert der geringsten Eindringtiefe o. Mittel
% minKonf=min(abs(N.k)); 
% rq=median(N.r(find(abs(N.k)==minKonf)));
rq=median(N.r);
if rq>20, rq=round(rq/10)*10; end
if ~iscell(Model.M), % grid model
    Model.M(:)=rq;
else
    for k=1:length(Model.M), Model.M{k}(:)=rq; end
end
Model.Bg=ones(size(Model.z))*rq;
Model.R=rq*ones(length(N.r),1);
Model.isfor=1;
messg(sprintf('Resetting Model to halfspace of %d Ohmm',round(rq)));
RMS=rms(N.r,Model.R,INV.lolo);
CHIQ=chi2(N.r,Model.R,N.err,INV.lolo);
dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
if ~isempty(datfile)
    [path,name,ext]=fileparts(datfile);
    if ~iscell(Model.M), % grid model
        matfile=strrep(datfile,ext,'-sens.mat');
    else % par model
        matfile=strrep(datfile,ext,'-s.mat');
    end
    if exist(matfile,'file'),
        S=[];
        load(matfile);
        global libmmfile
        libmmfile=checklic3;
        messg('Loading Sensitivity from disk');
        if iscell(Model.M), Coverage=zeros(Model.ncells,1); else Coverage=Model.M; end
        if issparse(S), % mirrored&sparse
            for i=1:size(S,1), Coverage(i)=sum(abs(S(i,:))); end
        else
            for i=1:size(S,2), Coverage(i)=sum(abs(S(:,i))); end
        end        
    else
        messg('Could not load Sensitivity!');
    end
end
if (INV.auto==0)&(INV.const==0), INV.lam=0; end %?????
if iscell(Model),
  Mref=[];for k=1:length(Model.M), Mref=[Mref;Model.M{k}(:)]; end
else Mref=Model.M(:); end
showback(handles,Model.Bg,Model.z);
dc3dinvres('showmod_Callback',gcbo,[],guidata(gcbo))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function start_Callback(h,eventdata,handles,varargin)
global Model N INV MAL FOR RMS CHIQ S
set(gcf,'Pointer','watch');
% [Model,dM]=onedinv(N,Model);
oldR=Model.R;dM=ones(size(S,2),1);
[Model.Bg,Model.R]=full1dinv(N,Model.z);
if iscell(Model.M),
    for k=1:length(Model.Bg)-1, Model.M{k}(:)=Model.Bg(k); end
else
    for k=1:length(Model.Bg)-1, Model.M(:,:,k)=Model.Bg(k); end
end
dc3dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
% set(handles.status,'String','Solving 1D forward problem...');
Model.R=fwd3d1d(N,Model.Bg,Model.z);
Model.isfor=1;
set(gcf,'Pointer','arrow');
set(handles.status,'String','Ready');
RMS=[RMS rms(N.r,Model.R,INV.lolo)];
CHIQ=[CHIQ chi2(N.r,Model.R,N.err,INV.lolo)];
dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))

% if ~isequal(size(S),[length(N.r),length(dM(:))]),
%     dc3dinvres('calcsens_Callback',gcbo,[],guidata(gcbo))
% end

% if isfield(INV,'broyden')&&(INV.broyden>0),
%     messg('Sensitivity updating by broyden method');
%     if isfield(INV,'lolo')&&(INV.lolo==0),
%         ddr=(Model.R-oldR-S*dM(:))/(dM(:)'*dM(:));ddr=ddr(:);
%     else
%         ddr=(log(Model.R)-log(oldR)-S*dM(:))/(dM(:)'*dM(:));ddr=ddr(:);
%     end
%     if issparse(S),
%        for i=1:size(S,2), 
%            fi=find(S(:,i));
%            S(fi,i)=S(fi,i)+ddr(fi)*dM(i); 
%        end
%     else
%       for i=1:size(S,2), S(:,i)=S(:,i)+ddr*dM(i); end
%     end
% end
if (INV.auto==0)&(INV.const==0), INV.lam=0; end % ???ex method
showback(handles,Model.Bg,Model.z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function step_Callback(h,eventdata,handles,varargin)
global S INV FOR N MAL RMS CHIQ Model Mref Coverage
if ~iscell(Model.M),
    Model.ncells=prod(size(Model.M));
else
    nn=0;for k=1:length(Model.M), nn=nn+prod(size(Model.M{k})); end
    Model.ncells=nn;
end
if ~isequal(sort(size(S)),sort([length(N.r),Model.ncells])),
    dc3dinvres('calcsens_Callback',gcbo,[],guidata(gcbo))
end

lbound=0;ubound=0;
if isfield(INV,'lbound'), lbound=INV.lbound; end
if isfield(INV,'ubound'), ubound=INV.ubound; end
if lbound>min(min(N.r),min(Model.R)), lbound=0; end
if (ubound>0)&&(ubound<max(max(N.r),max(Model.R))), ubound=0; end
set(handles.status,'String','Solving inverse subproblem...');

%% Calculating Residuum
if INV.lolo,
    dR=log(N.r(:)-lbound)-log(Model.R(:)-lbound);
    if ubound>0, dR=dR-log(ubound-N.r(:))+log(ubound-Model.R(:)); end
else
    dR=N.r(:)-Model.R(:);
end
% if INV.lolo==1, dR=log(N.r(:)-lbound)-log(Model.R(:)-lbound);
% else dR=N.r(:)-Model.R(:); end
% inverse subproblem
dBg=[]; % to fix!!!
if iscell(Model.M),
    MM=[];for k=1:length(Model.M), MM=[MM;Model.M{k}(:)]; end
else MM=Model.M(:); end
if INV.lolo,
    dM0=log(MM-lbound);
    if ubound>0, dM0=dM0-log(ubound-Model.M(:)); end
    if isequal(size(MM),size(Mref)),
        dM0=dM0-log(Mref-lbound);
        if ubound>0, dM0=dM0+log(ubound-Mref(:)); end
    end
else
    if isequal(numel(MM),numel(Mref)), dM0=MM(:)-Mref(:); end
end
set(gcf,'Pointer','watch');
[dM,lam]=minvers(S,dR,INV,N.err,Model,dM0);
tauopt=1;
oldR=Model.R;
if isfield(INV,'linesearch')&&(INV.linesearch>0),
    set(gcf,'Pointer','watch');
    newModel=modelupdate(Model,dM,2-INV.lolo,lbound,ubound);
    %if ~iscell(newModel.M), minmax(newModel.M), end
    set(figure(1),'MenuBar','none','NumberTitle','off','Name','Model');
    iconify(1);
    [cmin,cmax]=draw3dmodel(newModel,MAL,[],Coverage);drawnow;
    malcbar(handles,cmin,cmax,MAL.clog);
    Model.R=mfdfwd3d(newModel,N,FOR);
    [tauopt,appR]=linesearch(N,oldR,Model.R,INV.lolo);
    if (tauopt==0)||(chi2(N.r,Model.R,N.err,1)>CHIQ(end)*0.95), % line search failed
        taug=0.3; % sampling point
        newModel=modelupdate(Model,dM*taug,2-INV.lolo,lbound,ubound);
        Ri=mfdfwd3d(newModel,N,FOR);
        chiq=[CHIQ(end) chi2(N.r,Ri,N.err,1) chi2(N.r,Model.R,N.err,1)];
        taus=[0 taug 1];
        G=[1 0 0;1 taug taug^2;1 1 1];
        xx=G\chiq';
        tauopt=-xx(2)/2/xx(3);
%         taui=0:0.05:1; % plot results
%         chiqi=xx(1)+taui*xx(2)+taui.^2*xx(3);
%         plot(taui,chiqi,'bx-',taus,chiq,'ro',tauopt,xx(1)+tauopt*xx(2)+tauopt^2*xx(3),'g*');
        if tauopt<=0,    
            messg('(parabolic) line search failed, stopping iteration');
            dc3dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
            set(gcf,'Pointer','arrow');
            RMS(end+1)=RMS(end);
            CHIQ(end+1)=CHIQ(end);
            dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
            return;
        else
            if ~isfield(INV,'lolo')||(INV.lolo>0), 
                appR=oldR.*exp(tauopt*(log(Model.R)-log(oldR))); 
            else appR=oldR+tauopt*(Model.R-oldR); end
            messg(sprintf('Parabolic line search parameter %.2f, chi^2=%.1f-%.1f',...
                tauopt,chi2(N.r,Model.R,N.err,1),chi2(N.r,appR,N.err,1)));
        end
        if abs(tauopt-taug)<=0.5, % close to sample
            tauopt=taug; % take sampling point and save time
            dM=dM*tauopt;
            Model.R=Ri; 
            tauopt=1; % avoid new calculation
        end
    else
        if tauopt<1, messg(sprintf('Line search parameter %.2f, chi^2=%.1f-%.1f',...
                tauopt,chi2(N.r,Model.R,N.err,1),chi2(N.r,appR,N.err,1))); end
        dM=dM*tauopt;
    end
end
set(gcf,'Pointer','arrow');
% Model update
if INV.lolo==1, % logarithmic model parameters
    %     ma=find(dM==max(dM));ma=ma(1);mi=find(dM==min(dM));mi=mi(1);
    %     messg(sprintf('Max(dM) = %g Min(dM) = %g',exp(dM(ma)),exp(dM(mi))));
    Model=modelupdate(Model,dM,1,lbound,ubound);
    if length(dBg)==length(Model.Bg), Model.Bg(:)=Model.Bg(:).*exp(dBg(:)); end
else
    %     messg(sprintf('Max(dM) = %g Min(dM) = %g',max(dM),min(dM)));
    Model=modelupdate(Model,dM,2);
    if length(dBg)==length(Model.Bg), Model.Bg(:)=Model.Bg(:)+dBg(:); end
end
if iscell(Model.M),
    MM=[];for k=1:length(Model.M), MM=[MM;Model.M{k}(:)]; end
    mi=min(MM);ma=max(MM);
else
    mi=min(Model.M(:));ma=max(Model.M(:));
end
messg(sprintf('Min(Model) = %.1f Ohmm, Max(Model) = %.1f Ohmm',mi,ma));
%if ~iscell(Model.M), minmax(Model.M), end
dc3dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
drawnow;
set(handles.status,'String','Solving FD forward problem');
if FOR.method==2,
    rho=Model.Bg(1);
    if iscell(Model.M), % Para Model
        MM=[];for k=1:length(Model.M), MM=[MM;Model.M{k}(:)]; end
    else MM=Model.M(:); end
    Model.R=rho.*exp(S*(log(MM)-log(rho)));    
else
    set(gcf,'Pointer','watch');
    if ~isfield(INV,'linesearch')||(INV.linesearch==0)||(tauopt<0.95), Model.R=mfdfwd3d(Model,N,FOR); end
end
RMS=[RMS rms(N.r,Model.R,INV.lolo)];
CHIQ=[CHIQ chi2(N.r,Model.R,N.err,INV.lolo)];
dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
if isfield(INV,'robust')&&(INV.robust>0),
    if ~isfield(INV,'lolo')||(INV.lolo==0),
        dR=(N.r-Model.R)./N.err;
    else
        dR=(log(N.r)-log(Model.R))./log(1+N.err);
    end
    w=abs(dR)*sum(abs(dR))/sum(dR.^2);
    w(w<1)=1;
    N.err=N.err.*w;
    messg('Iteratively data reweighting (robust inversion)');
end
if ~isfield(INV,'broyden')||(INV.broyden>0),
    messg('Sensitivity updating by broyden method');
    if size(S,1)==numel(dM),
        if isfield(INV,'lolo')&&(INV.lolo==0),
            ddr=(Model.R-oldR-(dM(:)*S)')'/(dM(:)'*dM(:));
        else
            ddr=(log(Model.R)-log(oldR)-(dM(:)'*S)')'/(dM(:)'*dM(:));
        end
        if issparse(S),
            for i=1:size(S,1),
                fi=find(S(i,:));
                S(i,fi)=S(i,fi)+ddr(fi)*dM(i);
            end
        else
            for i=1:size(S,2), S(i,:)=S(i,:)+ddr*dM(i); end
        end
    else
        if isfield(INV,'lolo')&&(INV.lolo==0),
            ddr=(Model.R-oldR-S*dM(:))/(dM(:)'*dM(:));ddr=ddr(:);
        else
            ddr=(log(Model.R)-log(oldR)-S*dM(:))/(dM(:)'*dM(:));ddr=ddr(:);
        end
        if issparse(S),
            for i=1:size(S,2),
                fi=find(S(:,i));
                S(fi,i)=S(fi,i)+ddr(fi)*dM(i);
            end
        else
            for i=1:size(S,2), S(:,i)=S(:,i)+ddr*dM(i); end
        end
    end
end
Model.isfor=1;
set(handles.status,'String','Ready');
set(gcf,'Pointer','arrow');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function inversion_Callback(h,eventdata,handles,varargin)
global RMS Model INV MAL CHIQ
%invauto=INV.auto;INV.auto=1;
good=1;
while good,
    oldMod=Model;oldR=Model.R;
    dc3dinvres('step_Callback',gcbo,[],guidata(gcbo))
    dchiq=1-CHIQ(end)/CHIQ(end-1);
    if dchiq<0, % new model not better
        good=0;
        messg('Going back to old model (increasing CHI^2) and stop');
        Model=oldMod;Model.R=oldR;Model.isfor=1;
    end
    if CHIQ(end)<1, messg('CHI^2<1, stopping inversion'); end
    good=good&(CHIQ(end)>1)&(dchiq>0.05);
end
%INV.auto=invauto;

dc3dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
drawnow;

% --------------------------------------------------------------------
function calcsens_Callback(h, eventdata, handles, varargin)
global Model N S datfile INV Coverage
set(handles.status,'String','Building Sensitivity...');
[path,name,ext]=fileparts(datfile);
set(gcf,'Pointer','watch');
if iscell(Model.M), % para model
    S=modsens(Model,N);
    matfile=strrep(datfile,ext,'-s.mat');
    if ispc&&exist(matfile,'file'), dos(['attrib -H +A "' matfile '"']); end
    zsave=Model.z;
    save(matfile,'S','zsave');
else
    if isfield(INV,'spsens')&&(INV.spsens>0),
        t0=clock;
        S=spcalcsens3dt(Model.x,Model.y,Model.z,N,INV.spsens);
        messg(sprintf('ready(%ds) Elements: %s/%s(%.1f%%), %sB\n',...
            round(etime(clock,t0)),int2hum(nnz(S)),int2hum(prod(size(S))),...
            nnz(S)/prod(size(S))*100));
    else
        t0=clock;
        S=calcsens3d(Model.x,Model.y,Model.z,N);
        messg(sprintf('ready(%ds) %s elements (%sB)',...
            round(etime(clock,t0)),int2hum(prod(size(S))),int2hum(round(prod(size(S))*8))));
    end
    matfile=strrep(datfile,ext,'-sens.mat');
    if ispc&&exist(matfile,'file'), dos(['attrib -H +A "' matfile '"']); end
    zsave=Model.z;xsave=Model.x;ysave=Model.y;
    save(matfile,'S','xsave','ysave','zsave');
end
if iscell(Model.M), Coverage=zeros(Model.ncells,1); else Coverage=Model.M; end
% if isfield(INV,'spsens')&&(INV.spsens>0), % mirrored&sparse
if issparse(S), % mirrored&sparse
    for i=1:size(S,1), Coverage(i)=sum(abs(S(i,:))); end
else
    for i=1:size(S,2), Coverage(i)=sum(abs(S(:,i))); end
end
if ispc, dos(['attrib +H -A "' matfile '"']); end
set(gcf,'Pointer','arrow');
set(handles.status,'String','Ready');

% --------------------------------------------------------------------
function varargout = delsens_Callback(h, eventdata, handles, varargin)
global S
S=[];

% --------------------------------------------------------------------
function varargout = showlcurve_Callback(h, eventdata, handles, varargin)
load('rhoeta.mat');
if exist('rho','var')&exist('eta','var')&exist('ak','var'),
    set(figure(3),'MenuBar','none','NumberTitle','off','Name','L-Curve');
    iconify(fig);
    lkurv(rho,eta,ak);
end

% -------- FILE ------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = readdata_Callback(h, eventdata, handles, varargin)
global N S Model MAL RMS CHIQ datfile INV Coverage
alltypes='*.dat;*.pro;*.res;*.s3d;*.txt;*.gem;*.ohm';
if isempty(datfile), infile=alltypes; else
    [pname,fname,ext]=fileparts(datfile);
    infile=strrep(datfile,[fname ext],alltypes);
end 
[filename,pathname]=uigetfile({alltypes,'Known Types';'*.dat;*.ohm',...
    'DAT-files (*.dat/ohm)';'*.pro','profiles files (*.pro)';...
    '*.s3d','3d soundings files (*.s3d)';...
    '*.gem','sensinv3d files (*.gem)';...
    '*.txt','resecs files (*.txt)';...
    '*.*','All Files (*.*)'},'Read datafile',infile);
if isequal(filename,0), return; end
datfile=strrep(fullfile(pathname,filename),[pwd filesep],'');
[fpath,name,ext]=fileparts(datfile);
global libmmfile
libmmfile=checklic3;
if isequal(lower(ext),'.s3d'),
    N=readsnd3d(datfile);
elseif isequal(lower(ext),'.pro')
    N=readpro(datfile);
%     if any(N.elec(:,3)~=0),
%         N.topo=N.elec;
%         N.elec(:,3)=0;
%     end
elseif isequal(lower(ext),'.gem'),
    N=readgemfile(datfile);
elseif isequal(lower(ext),'.txt'),
    N=readresecsfile(datfile);
else
    N=read3dfile(datfile);
end
delmeas=[];
if isfield(N,'u'), delmeas=find(N.u==0); end
if isfield(N,'rho'), delmeas=[delmeas;find(N.rho==0)]; end
if isfield(N,'r'), delmeas=[delmeas;find(N.r==0)]; end
if ~isempty(delmeas),
    delmeas=unique(delmeas);
    uiwait(errordlg(sprintf('Found %d zero values. Deleting.',length(delmeas)),'Zero values'));
    N=delmeasurement(N,delmeas);
end
if (~isfield(N,'elec'))||(~isfield(N,'a')),
    messg('Improper data structure');return; end
if isfield(N,'d'),
    N.elec(:,3)=abs(N.d);
    N.k=getkonf3d(N);    
else
    if (size(N.elec,2)>2)&&any(N.elec(:,3)), % topo found
        N.topo=N.elec; %save topo
        expfirst=questdlg('Neglecting topography! Export ohm-File before?');
        if isequal(expfirst,'Yes'),
            if ~isfield(N,'rho'), N.rho=N.r./N.k; end
            [fname,pname]=uiputfile(strrep(datfile,ext,'.ohm'),'Export Ohm-File');
            if ischar(fname), 
                outfile=fullfile(pname,fname);
                if ~isfield(N,'rho')||isempty(N.rho), N.rho=N.r./N.k; end
                saveinv3dfile(outfile,rmfield(N,'r')); 
            end
            N.elec(:,3)=0;
        elseif isequal(expfirst,'No'),
            N.elec(:,3)=0;
        else
            if isfield(N,'zweid'), N=rmfield(N,'zweid'); end
        end
        if isfield(N,'zweid'),
            for i=1:length(N.zweid),
                el=N.zweid{i}.elec;
                el(2:end,1)=cumsum(round(sqrt(sum(diff(el).^2,2))*10)/10);
                el(1,1)=0;el(:,2)=0;
                N.zweid{i}.elec=el;
            end
        end
    end
end

un=unique(N.b);
if (length(un)==1)&&(un>0),
    if isequal(questdlg('Treat constant current electrode as infinite?'),'Yes'),
        N.b(:)=0;
%         ff=1:size(N.elec,1);ff(un:end)=ff(un:end)+1;N.a=ff(N.a);N.b=ff(N.b);
    end 
end
un=unique(N.n);
if (length(un)==1)&&(un>0),
    if isequal(questdlg('Treat constant potential electrode as infinite?'),'Yes'),
        N.n(:)=0; end 
end
if ~isfield(N,'k'), N.k=getkonf3d(N); end
if ~isfield(N,'r'),
    if isfield(N,'rho'), N.r=N.rho.*N.k;
    elseif isfield(N,'u')&&isfield(N,'i'), N.r=N.u./N.i.*N.k;
    else errordlg('Could not derive apparent resistivity');return; end
end    
[fi1,fi2]=getrez(N);
if ~isempty(fi1),
    messg(sprintf('Found reciprocals for %d data',length(fi1)));
    R1=N.r(fi1);
    R2=N.r(fi2);
    set(figure(9),'MenuBar','none','NumberTitle','off','Name','Reciprocity crossplot');
    plot(R1,R2,'.');xlabel('normal');ylabel('reverse');grid on;
%     mima=minmax(N.r);
    mima=[min(min(R1),min(R2)) max(max(R1),max(R2))];
    set(gca,'Xlim',mima,'Ylim',mima);line(mima,mima);
    rez=(R1-R2)./(R1+R2);
    messg(sprintf('std/max reciprocity = %.1f/%.1f%%',...
        std(rez)*100,max(rez)*100));
    N.rez=N.a*0;
    N.rez(fi1)=rez;
    N.rez(fi2)=rez;
    drawnow;
end
if isempty(N), uiwait(errordlg('File format unknown!'));return; end
if (~isfield(N,'err'))||(length(N.err)~=length(N.a))||(min(N.err)<=0),
    daterr;
else
    messg(sprintf('Found error estimate in file, min=%.1f%%, max=%.1f%%',min(N.err)*100,max(N.err)*100));
    if min(N.err)<0.01,
        aa=inputdlg('Type percentage error to be added (2%)','Small errors found!');
        bb=str2num(aa{1});
        if ~isempty(bb), bb=2; end
        N.err=N.err+bb/100;
    end
end
[pa,name,ext]=fileparts(datfile);
if (~isfield(N,'zweid'))&&(~isfield(N,'eind'))&&(~any(N.elec(:,3))),
    N=getpseudos(N); % try to get pseudosections
end
if length(N.r)<1000,
  dc3dinvres('showdata_Callback',gcbo,[],guidata(gcbo)); end
% Creating model
S=[];ss='Grid'; % grid model default
% ss=questdlg('Grid model or parametric model?','Choose Model Type','Grid','Para','Grid');
if strcmp(ss,'Para'),
    matfile=strrep(datfile,ext,'-s.mat');
    Model=create3dmod(N);
    if exist(matfile,'file'), 
        load(matfile); 
        global libmmfile
        libmmfile=checklic3;
        messg('Loading Sensitivity matrix from disk');
    end
else % grid model
    Model=[];
    matfile=strrep(datfile,ext,'-sens.mat');
    if exist(matfile,'file'), 
        xsave=[];ysave=[];zsave=[];
        load(matfile); 
        global libmmfile
        libmmfile=checklic3;
        messg('Loading Sensitivity matrix from disk');
        if any(xsave)&any(ysave)&any(zsave),
            Model.x=xsave;Model.y=ysave;Model.z=zsave;
            rq=median(N.r);if rq>10, rq=round(rq); end
            Model.M=ones(length(Model.x)-1,length(Model.y)-1,length(Model.z)-1)*rq;
            Model.Bg=ones(size(Model.z))*rq;
            messg(sprintf('Loading %dx%dx%d=%d cell model from sensitivity',...
                size(Model.M,1),size(Model.M,2),size(Model.M,3),prod(size(Model.M))));
        else
            %[Model.x,Model.y,Model.z,Model.M,Model.Bg]=modelfromdata3d(N);    
            Model=modelfromdata3d(N);    
        end
    else % ~ exist matfile
        Model=modelfromdata3d(N);    
    end
    Model.ncells=prod(size(Model.M));
    if iscell(Model.M), 
        Coverage=zeros(Model.ncells,1); 
    else
        Coverage=zeros(size(Model.M));
    end
    if issparse(S), % mirrored&sparse
        for i=1:size(S,1), Coverage(i)=sum(abs(S(i,:))); end
    else
        for i=1:size(S,2), Coverage(i)=sum(abs(S(:,i))); end
    end        
end
if isempty(S),
    space=round(Model.ncells*length(N.r)*8/1024/1024);
    if space>800, %Megabytes
        warndlg(['Attention! The model will cause a sensitivity matrix of ' num2str(space)...
                'MB RAM. Consider changing the model or using sparse sensitivity!'],'Model too large!');
    else
        messg(sprintf('Estimated size of (full) sensitivity matrix=%dMB',space));
    end
end
Model.R=Model.Bg(1)*ones(length(N.r),1);
Model.isfor=1;

if (INV.auto==0)&(INV.const==0), INV.lam=0; end %?? ex method
showback(handles,Model.Bg,Model.z);
RMS=rms(N.r,Model.R,INV.lolo);
CHIQ=chi2(N.r,Model.R,N.err,INV.lolo);
dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
dc3dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
dc3dinvres('clearsvd_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function varargout = readpro_Callback(h, eventdata, handles, varargin)
global N S Model MAL RMS CHIQ datfile INV Coverage
% Reading & plotting data, estimating error
if isempty(datfile), infile='*.pro'; else
    [pname,fname,ext]=fileparts(datfile);
    infile=strrep(datfile,[fname ext],'*.pro');
end 
[filename,pathname]=uigetfile({'*.pro','PRO-files (*.pro)';'*.*','All Files (*.*)'},'Read PRO File',infile);
if isequal(filename,0), return; end
ss=get(handles.output,'String');
sss={};
for l=1:4, sss{l}=ss{l}; end
set(handles.output,'String',sss);
datfile=strrep(fullfile(pathname,filename),[pwd filesep],'');
N=readpro(datfile);
if isempty(N)||~isfield(N,'r')||(length(N.r)==0), 
    uiwait(errordlg('PRO-File incorrupt!','Bad pro file'));
    return; 
end
if (size(N.elec,2)>2)&&(min(N.elec(:,3))~=max(N.elec(:,3))),
    N.elec(:,3)=0;
    messg('Found topography in data files, setting all electrodes to zero');
end
dc3dinvres('showdata_Callback',gcbo,[],guidata(gcbo))
if (~isfield(N,'err'))||(length(N.err)==0)||(min(N.err)<=0), daterr; end
[path,name,ext]=fileparts(datfile);
S=[];
% Creating & plotting Model
if strcmp(questdlg('Grid model or parametric model?',...
        'Choose Model Type','Grid','Para','Grid'),'Grid'),
    Model=[];
    matfile=strrep(datfile,ext,'-sens.mat');
    if exist(matfile,'file'), 
        global libmmfile
        libmmfile=checklic3;
        load(matfile); 
        messg('Loading Sensitivity matrix from disk');
        if exist('xsave','var')&exist('ysave','var')&exist('zsave','var'),
            Model.x=xsave;Model.y=ysave;Model.z=zsave;
            rq=median(N.r);if rq>10, rq=round(rq); end
            Model.M=ones(length(Model.x)-1,length(Model.y)-1,length(Model.z)-1)*rq;
            Model.Bg=ones(size(Model.z))*rq;
            messg(sprintf('Loading %dx%dx%d=%d cell model from sensitivity',...
                size(Model.M,1),size(Model.M,2),size(Model.M,3),prod(size(Model.M))));
        else
            [Model.x,Model.y,Model.z,Model.M,Model.Bg]=modelfromdata3d(N);    
        end
    else
        [Model.x,Model.y,Model.z,Model.M,Model.Bg]=modelfromdata3d(N);    
    end
    Model.ncells=prod(size(Model.M));
    if iscell(Model.M), Coverage=zeros(Model.ncells,1); else Coverage=Model.M; end
    if issparse(S), % mirrored&sparse
        for i=1:size(S,1), Coverage(i)=sum(abs(S(i,:))); end
    else
        for i=1:size(S,2), Coverage(i)=sum(abs(S(:,i))); end
    end        
else %para model
    matfile=strrep(datfile,ext,'-s.mat');
    Model=create3dmod(N);
    if exist(matfile,'file'), 
        global libmmfile
        libmmfile=checklic3;
        load(matfile); 
        messg('Loading Sensitivity matrix from disk');
    end
end
Model.R=Model.Bg(1)*ones(length(N.r),1);
Model.isfor=1;
showback(handles,Model.Bg,Model.z);
RMS=rms(N.r,Model.R,INV.lolo);
CHIQ=chi2(N.r,Model.R,N.err,INV.lolo);
dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
dc3dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
dc3dinvres('clearsvd_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function load_Callback(h,eventdata,handles,varargin)
global S N MAL INV dR Model RMS datfile CHIQ Coverage
[path,name,ext]=fileparts(datfile);
if isempty(datfile)
    matfile='*.mat';
else
    matfile=strrep(datfile,ext,'.mat');
    if ~exist(matfile,'file'), matfile='*.mat'; end
end
[fname,pname]=uigetfile(matfile,'Load Workspace');
if isequal(fname,0), return; end
matfile=fullfile(pname,fname);
Model=[];
ohandles=handles;
global libmmfile
libmmfile=checklic3;
load(matfile);
handles=ohandles;
if isempty(Model), % old grid model
    if exist('x','var'), Model.x=x; end
    if exist('y','var'), Model.y=y; end
    if exist('z','var'), Model.z=z; end
    if exist('Bg','var'), Model.Bg=Bg; end
    if exist('M','var'), 
        Model.M=M; 
        Model.ncells=prod(size(M));
    end
end
if ~isfield(N,'err'), N.err=estimateerror(N); end
if isempty(CHIQ), CHIQ=chi2(N.r,Model.R,N.err,INV.lolo); end
if ~isfield(N,'zweid'), N=getpseudos(N); end
if exist('mess','var')&&(~isempty(mess))&ishandle(handles.output), 
    set(handles.output,'String',mess); 
else
    messg('(Loading workspace from disk file)');
    dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
end
%set(gcbo,'name',[ 'INV2D - ' datfile ]);
[path,name,ext]=fileparts(datfile);
if iscell(Model.M),
    matfile=strrep(datfile,ext,'-s.mat');
else
    matfile=strrep(datfile,ext,'-sens.mat');
end
if isempty(S)&(exist(matfile,'file')), 
    load(matfile);
    global libmmfile
    libmmfile=checklic3;
    messg('Loaded Sensitivity matrix from disk');
    if iscell(Model.M), Coverage=zeros(Model.ncells,1); else Coverage=Model.M; end
    if issparse(S), % mirrored&sparse
        for i=1:size(S,1), Coverage(i)=sum(abs(S(i,:))); end
    else
        for i=1:size(S,2), Coverage(i)=sum(abs(S(:,i))); end
    end        
    if (~isequal(z,zsave)),
        S=[];
        messg('Recomputation of sensitivity will be necessary');
    end
end
dc3dinvres('showmod_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------
function save_Callback(h,eventdata,handles,varargin)
global S N MAL INV FOR Model RMS datfile CHIQ
mess=get(handles.output,'String');
if isempty(datfile), matfile='*.mat'; else
    [path,name,ext]=fileparts(datfile);
    matfile=strrep(datfile,ext,'.mat'); end
[fname,pname]=uiputfile(matfile,'Save Workspace');
if fname~=0, % Testversion n. Zeile einklammern
    global libmmfile
    if ~isequal(libmmfile,4), errordlg('Save workspace not supported in Test version!');return; end
    matfile=fullfile(pname,fname);
    clear path name ext fname pname
    messg(['Saving workspace to file ',matfile]);
    save(matfile);
end

% --------------------------------------------------------------------
function varargout = adddata_Callback(h, eventdata, handles, varargin)
global Model N RMS CHIQ datfile INV
alltypes='*.dat;*.pro;*.res;*.s3d;*.txt;*.gem;*.ohm';
if isempty(datfile), infile=alltypes; else
    [pname,fname,ext]=fileparts(datfile);
    infile=strrep(datfile,[fname ext],'*.dat');
end 
[filename,pathname]=uigetfile({alltypes,'Known Types';'*.dat;*.ohm',...
    'DAT-files (*.dat/ohm)';'*.pro','profiles files (*.pro)';...
    '*.s3d','3d soundings files (*.s3d)';...
    '*.gem','sensinv3d files (*.gem)';...
    '*.txt','resecs files (*.txt)';...
    '*.*','All Files (*.*)'},'Read datafile',infile);
if isequal(filename,0), return; end
datfile=strrep(fullfile(pathname,filename),[pwd filesep],'');
N1=N;
if isequal(lower(ext),'.s3d'),
    N=readsnd3d(datfile);
elseif isequal(lower(ext),'.pro')
    N=readpro(datfile);
elseif isequal(lower(ext),'.txt'),
    N=readresecsfile(datfile);
else
    N=read3dfile(datfile);
end
% switch check3dfile(datfile),
% case 1,
%     N=readres3dinvfile(datfile);
% case 2,
%     N=readinv3dfile(datfile);
% case 3,
%     N=readres3dinvfile(datfile,1);
% otherwise
%     messg('COuld not determine data type');
%     return;
% end
if (~isfield(N,'err'))||isempty(N.err)||(min(N.err)<0),
    daterr;
else
    messg(sprintf('Found error estimate in file, min=%.1f%%, max=%.1f%%',min(N.err)*100,max(N.err)*100));
    if min(N.err)<0.01,
        aa=inputdlg('Type percentage error to be added (2%)','Small errors found!');
        bb=str2num(aa{1});
        if isempty(bb), bb=2; end
        N.err=N.err+bb/100;
    end
end
N=combdata3d(N1,N);
N=getpseudos(N);
if ~isfield(Model,'M')||isempty(Model.M), 
    Model=create3dmod(N); 
    dc3dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
end
Model.R=Model.Bg(1)*ones(length(N.r),1);
Model.isfor=1;
S=[];
% [path,name,ext]=fileparts(datfile);
% matfile=strrep(datfile,ext,'-s.mat');
% if exist(matfile)==2, 
%     load(matfile); 
%     messg('Loading Sensitivity matrix from disk');
% end

if (INV.auto==0)&(INV.const==0), INV.lam=0; end %?? ex method
showback(handles,Model.Bg,Model.z);
RMS=rms(N.r,Model.R,INV.lolo);
CHIQ=chi2(N.r,Model.R,N.err,INV.lolo);
dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
dc3dinvres('clearsvd_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function varargout = clearmessages_Callback(h, eventdata, handles, varargin)
set(handles.output,'String',...
    {'DC3DInvRes - DC 3D Inversion & Resolution','Thomas Günther - RESISTIVITY.NET',...
    'Version: 2.11.6','----------------------------------------------'});

% --------------------------------------------------------------------
function varargout = savemessages_Callback(h, eventdata, handles, varargin)
fid=fopen('dc3dinvres.log','w');
ss=get(handles.output,'String');
fprintf(fid,'%s\n',ss{:});
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ende_Callback(h,eventdata,handles,varargin)
global output datfile FOR INV MAL
output=[];ff=datfile;
if exist('setup.mat','file'), load('setup.mat'); end
datfile=ff;
save('setup.mat','FOR','INV','MAL','datfile'); 
mm=[tempdir 'model.mat'];
if exist(mm,'file'), delete(mm); end
for i=1:10,
    if ishandle(i), close(i); end
end
delete(gcbf);


% ------- SHOW -------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showmod_Callback(h,eventdata,handles,varargin)
global Model MAL Coverage
set(figure(1),'MenuBar','none','NumberTitle','off','Name','Model');
iconify(1);
[cmin,cmax]=draw3dmodel(Model,MAL,[],Coverage);
malcbar(handles,cmin,cmax,MAL.clog);
showback(handles,Model.Bg,Model.z);
if isfield(handles,'status'), set(handles.status,'String','Ready'); end
if isfield(handles,'slidesens'), set(handles.slidesens,'Visible','Off'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = showcov_Callback(h, eventdata, handles, varargin)
global MAL S Model Coverage
set(figure(1),'MenuBar','none','NumberTitle','off','Name','Coverage');
iconify(1);
mal=MAL;
mal.clog=0;
mal.cauto=1;
mal.alpha=0;
[cmin,cmax]=draw3dmodel(modelupdate(Model,Coverage,0),mal);
malcbar(handles,cmin,cmax,0,'');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = showvol_Callback(h, eventdata, handles, varargin)
global x y z M Model
if iscell(Model.M), % para model
    [M,x,y,z]=mesch3dmodel(Model);
else
    x=Model.x;y=Model.y;z=Model.z;M=Model.M;
end
vol3d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showback(handles,Bg,Z)
if isequal(Bg,0),
    old='';
else
    old=sprintf('Background\nresistivities:\nz/m rho/Om\n%4.1f\n');
    for l = 1:min(length(Bg),length(Z)),
        old=sprintf('%s%4.1f:\n%10d\n',old,Z(l),round(Bg(l)));
    end
end
if isfield(handles,'bg')&&ishandle(handles.bg), set(handles.bg,'String',old); end

% --------------------------------------------------------------------
function varargout = showdata_Callback(h, eventdata, handles, varargin)
global N MAL
mal=MAL;mal.vie=1;
set(figure(2),'MenuBar','none','NumberTitle','off','Name','Data');
iconify(2);
if mal.cauto,
    mal.cmin=min(N.r);mal.cmax=max(N.r);
    mal.cauto=0;mal.clog=(mal.cmin>0);
    aa=sort(N.r);mal.cmax=aa(fix(length(aa)*0.95));
end
if isfield(N,'zweid')&&~isempty(N.zweid),
    [cmin,cmax]=plotprofiles(N,[],mal);
    malcbar(handles,cmin,cmax,mal.clog,'\rho_a in \Omega\cdotm');
else
    loghist(N.r,100);
end

% --------------------------------------------------------------------
function varargout = showresponse_Callback(h, eventdata, handles, varargin)
global N Model MAL
set(figure(3),'MenuBar','none','NumberTitle','off','Name','Response');
iconify(3);
mal=MAL;mal.vie=1;
if mal.cauto,
    mal.cmin=min(N.r);mal.cmax=max(N.r);
    aa=sort(N.r);mal.cmax=aa(fix(length(aa)*0.95));
    mal.cauto=0;mal.clog=(mal.cmin>0);
end
[cmin,cmax]=plotprofiles(N,Model.R,mal);
malcbar(handles,cmin,cmax,mal.clog,'\rho_a in \Omega\cdotm');

% --------------------------------------------------------------------
function varargout = showmisfit_Callback(h, eventdata, handles, varargin)
global N Model MAL
set(figure(2),'MenuBar','none','NumberTitle','off','Name','Misfit');
iconify(2);
mal=MAL;mal.vie=1;mal.clog=0;mal.clog=0;mal.cauto=1;mal.cmap=2;
dR=(1-Model.R./N.r)*100;
%mm=max(abs(dR));mal.cmin=-mm;mal.cmax=mm;
[cmin,cmax]=plotprofiles(N,dR,mal);
malcbar(handles,cmin,cmax,mal.clog,'\Delta in %');

% --------------------------------------------------------------------
function showerror_Callback(hObject, eventdata, handles)
global N Model MAL
set(figure(2),'MenuBar','none','NumberTitle','off','Name','Error');
iconify(2);
mal=MAL;mal.vie=1;mal.clog=1;mal.cauto=1;
[cmin,cmax]=plotprofiles(N,N.err*100,mal);
malcbar(handles,cmin,cmax,1,'\epsilon in %');

% ---------------------------------------
function malcbar(handles,cmin,cmax,clog,titel)
if nargin<5, titel='\rho in \Omega\cdotm'; end
if nargin<4, clog=0; end
% n=64;image([0 1],[1 n 1],(1:n)','parent',handles.cbar);
% cbar(cmin,cmax,clog,1);
% return;
if clog, cmin=log10(cmin);cmax=log10(cmax); end
cmap=colormap;
if isequal(cmap(1,:),[1 1 1]), cmap(1,:)=[]; end
n=size(cmap,1);
if ~isfield(handles,'cbar'), return; end
image([0 1],[1 n 1],(1:n)','parent',handles.cbar);
set(handles.cbar,'DataAspectRatio',[1 1 1]);
colormap(handles.cbar,cmap);
set(handles.cbar,'XTick',[],'YDir','normal');
set(handles.cbar,'YTick',linspace(1,n,11));
yt=get(handles.cbar,'YTick');
ytl=yt/n*(cmax-cmin)+cmin;
if clog, ytl=10.^ytl; end

% fi=find(ytl>=10);ytl(fi)=round(ytl(fi));
% fi=find(ytl<10);
% for i=1:length(fi),
%     l=0;yy=ytl(fi(i));
%     if yy~=0,
%         while abs(yy)<10, yy=yy*10;l=l+1; end
%         yy=round(yy);
%         for ll=1:l, yy=yy/10; end
%     end
%     ytl(fi(i))=yy;
% end
set(handles.cbar,'YTickLabel',num2strcell(rndig(ytl)));
set(get(handles.cbar,'Title'),'String',titel);

% ----------------- OPTIONS ------------------------------------------
% --------------------------------------------------------------------
function varargout = optmod_Callback(h, eventdata, handles, varargin)
global Model N FOR MAL RMS CHIQ INV S
oldBg=Model.Bg;
if iscell(Model.M), %para model
    m_mod;
    nn=0;
    for k=1:length(Model.M), nn=nn+prod(size(Model.M{k})); end
    Model.ncells=nn;
else
    ncells=Model.ncells;oldz=Model.z;oldx=Model.x;
    s_mod;
    if ~isequal(Model.z,oldz)|~isequal(Model.x,oldx), 
        S=[];         
        messg(sprintf('I=%d J=%d K=%d => %d Cells',length(Model.x)-1,length(Model.y)-1,length(Model.z)-1,numel(Model.M)));
        messg(sprintf('X=%.1f-%.1f Y=%.1f-%.1f Z=%.1f-%.1f',Model.x(1),Model.x(end),Model.y(1),Model.y(end),Model.z(1),Model.z(end)));
    end
    Model.ncells=prod(size(Model.M));
    if Model.ncells~=ncells,
        space=round(Model.ncells*length(N.r)*8/1024/1024);
        if space>800, %Megabytes
            warndlg(['Attention! The model will cause a sensitivity matrix of ' num2str(space)...
                'MB RAM. Consider changing the model or using sparse sensitivity!'],'Model too large!');
        else
            messg(sprintf('Estimated size of (full) sensitivity matrix=%dMB',space));
        end
    end
end
if ~isequal(oldBg,Model.Bg)&(length(unique(Model.Bg))>1),
    messg('Changed Background, calculating 1D forward response');
    set(handles.status,'String','Solving 1D forward problem...');
    set(gcf,'Pointer','watch');
    if length(unique(Model.Bg))>1,
        Model.R=fwd3d1d(N,Model.Bg,Model.z);
    else
        Model.R=ones(length(N.r),1)*Model.Bg(1);
    end
    Model.isfor=1;
    set(gcf,'Pointer','arrow');
    %Model.R=mfdfwd3d(Model,N,FOR);
    RMS=[RMS(1) rms(N.r,Model.R,INV.lolo)];
    CHIQ=[CHIQ(1) chi2(N.r,Model.R,N.err,INV.lolo)];
    dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
    set(handles.status,'String','Ready');
else
    if length(unique(Model.Bg))==1,
        Model.R=ones(size(Model.R))*Model.Bg(1);
        Model.isfor=1;
        RMS=rms(N.r,Model.R,INV.lolo);
        CHIQ=chi2(N.r,Model.R,N.err,INV.lolo);
        dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo));
    end
end
dc3dinvres('showmod_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function varargout = optinv_Callback(h, eventdata, handles, varargin)
s_inv;

% --------------------------------------------------------------------
function varargout = optfor_Callback(h, eventdata, handles, varargin)
m_for;

% --------------------------------------------------------------------
function varargout = optmal_Callback(h, eventdata, handles, varargin)
m_mal;

% --------------------------------------------------------------------
function varargout = optsave_Callback(h, eventdata, handles, varargin)
global INV FOR MAL
save('setup.mat','INV','MAL','FOR');

% --------------------------------------------------------------------
function varargout = optload_Callback(h, eventdata, handles, varargin)
global INV FOR MAL
if exist('setup.mat','file'),
    load('setup.mat');
else
    messg('Could not load Setup File!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function standard_Callback(h,eventdata,handles,varargin)
global INV MAL FOR
INV=struct('method',0,'redu',0,'mitschicht',0,'lolo',1,'lam',30,'mico',0.4,'sens',0,'auto',2,'update',2,'const',0,'weight',1,'lbound',0,'glob',1,'linesearch',1,'spsens',0);
FOR=struct('method',0,'acc',1e-4,'tol',1e-4,'maxit',50,'rand',4,'prolong',4,'zusatz',2,'refine',1,'fillup',1,'direct',-1);
MAL=struct('xy',1,'clog',1,'cauto',1,'cmap',5,'nu',0,'nv',0,'startwith',0,'elec',1,'xdir',0,'ydir',0,'vie',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function muchdata_Callback(h,eventdata,handles,varargin)
global INV MAL FOR
INV=struct('method',2,'redu',0,'mitschicht',0,'lolo',1,'lam',0.1,'mico',0.4,'sens',0,'auto',0);
FOR=struct('method',0,'acc',1e-3,'tol',1e-4,'maxit',50,'rand',3,'prolong',5,'zusatz',1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fastsirt_Callback(h,eventdata,handles,varargin)
global INV MAL FOR
INV=struct('method',4,'redu',0,'mitschicht',0,'lolo',1,'lam',0.1,'mico',0.4,'sens',0,'auto',1);
FOR=struct('method',2,'acc',1e-3,'tol',1e-4,'maxit',50,'rand',3,'prolong',5,'zusatz',1);

% ----------------- MODEL ------------------------------------------
% --------------------------------------------------------------------
function varargout = modfree_Callback(h, eventdata, handles, varargin)
modelling;

% --------------------------------------------------------------------
function varargout = modsave_Callback(h, eventdata, handles, varargin)
global Model
save([tempdir 'model.mat'],'Model');

% --------------------------------------------------------------------
function varargout = modload_Callback(h, eventdata, handles, varargin)
global Model
modfile=[tempdir 'model.mat'];
if exist(modfile,'file'),
    load(modfile);
else
    messg('Could NOT read model file!');
end

% --------------------------------------------------------------------
function varargout = modexp_Callback(h, eventdata, handles, varargin)
global Model datfile Coverage
[path,name,ext]=fileparts(datfile);
outfile=strrep(datfile,ext,'.mod');
[fname,pname]=uiputfile(outfile,'Export model to ASCII-File');
if fname==0, return; end% no valid name
outfile=fullfile(pname,fname);
[ff,pp,ee]=fileparts(outfile);
if strcmp(ee,''), outfile=[outfile '.mod']; end
messg(['Exporting model to ASCII File : ' outfile]);
if iscell(Model.M), % para model
    messg('Not yet possible for para model');
else % grid model
    modelexport3d(outfile,Model.M,Model.x,Model.y,Model.z,Coverage);
%     warndlg('Resistivity values are replaced by zero in test version!','Warning');
end % Testversion .M*0

% --------------------------------------------------------------------
function varargout = modimp_Callback(h, eventdata, handles, varargin)
global Model datfile Coverage N RMS CHIQ INV
outfile='*.mod';
if ~isempty(datfile),
    [path,name,ext]=fileparts(datfile);
    outfile=strrep(datfile,ext,'.mod');
    if ~exist(outfile,'file'), outfile=strrep(outfile,name,'*'); end
end
[fname,pname]=uigetfile({'*.MOD','model files (*.mod)';'*.*','All Files (*.*)'},'Read datafile',outfile);
%[fname,pname]=uigetfile(outfile,'Import model from ASCII-File');
if isequal(fname,0), return; end
outfile=fullfile(pname,fname);
if 1,
    [Model.M,Model.x,Model.y,Model.z,Covi]=modelimport3d(outfile);
    if max(Covi(:))>0, Coverage=Covi; end
    islay=1;
    for i=1:size(Model.M,3),
        mm=Model.M(:,:,i);
        Model.Bg(i)=median(mm(:));
        islay=islay&&(length(unique(mm(:)))==1);
    end
    Model.Bg(length(Model.z))=Model.Bg(size(Model.M,3));
    if islay, Model.R=fwd3d1d(N,Model.Bg,Model.z);
        Model.isfor=1; end
    RMS=rms(N.r,Model.R,INV.lolo);
    CHIQ=chi2(N.r,Model.R,N.err,INV.lolo);
    dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
%     Model.Bg=0;
else
    fid=fopen(outfile,'r'); 
    zeile='';
    while length(zeile)==0, zeile=destrip(fgetl(fid)); end
    XYZR=fscanf(fid,'%f',[7,Inf])';
    XYZR(end+1,:)=sscanf(zeile,'%f',[7,1])';
    x=unique(XYZR(:,1));x(end+1)=max(XYZR(:,2));
    y=unique(XYZR(:,3));y(end+1)=max(XYZR(:,4));
    z=unique(XYZR(:,5));z(end+1)=max(XYZR(:,6));
    Model.M=ones(length(x)-1,length(y)-1,length(z)-1);
    [tf,ii]=ismember(XYZR(:,1),x);
    [tf,jj]=ismember(XYZR(:,3),y);
    [tf,kk]=ismember(XYZR(:,5),z);
    for l=1:length(ii),
        Model.M(ii(l),jj(l),kk(l))=XYZR(l,7);
    end
    % while ischar(zeile),
    %     zeile='';
    %     while length(zeile)==0, zeile=destrip(fgetl(fid)); end
    % end
    fclose(fid);
    Model.x=x;Model.y=y;Model.z=z;Model.Bg=0;
end
messg(['Imported model from ASCII File : ' outfile]);
if isempty(datfile)||~isempty(strfind(datfile,'.mod')), datfile=outfile; end
dc3dinvres('showmod_Callback',gcbo,[],guidata(gcbo))

function zeile=destrip(zeile)
% strip string from comments (with # character)
aa=strfind(zeile,'#');
if ~isempty(aa), zeile=zeile(1:aa(1)-1); end

% --------------------------------------------------------------------
function varargout = modsavef_Callback(h, eventdata, handles, varargin)
global Model N
fid=fopen('modelling.out','w');
for l=1:length(N.r)
    aa=N.elec(N.a(l),1:2);
    fprintf(fid,'%.1f ',aa);
    if N.b(l)>0,
        aa=N.elec(N.b(l),1:2);
        fprintf(fid,'%.1f ',aa);
    end
    aa=N.elec(N.m(l),1:2);
    fprintf(fid,'%.1f ',aa);
    if N.n(l)>0,
        aa=N.elec(N.n(l),1:2);
        fprintf(fid,'%.1f ',aa);
    end
    fprintf(fid,'%.2f\n',Model.R(l));
end
fclose(fid);

% --------------------------------------------------------------------
function varargout = modopen_Callback(h, eventdata, handles, varargin)
global Model
[filename,pathname]=uigetfile({'*.3d'},'Read 3d-file');
datfile=fullfile(pathname,filename);

% --------------------------------------------------------------------
function varargout = output_Callback(h, eventdata, handles, varargin)


% --------------------------------------------------------------------
function varargout = figure1_KeyPressFcn(h, eventdata, handles, varargin)
switch get(h,'CurrentCharacter')
    case '1'
        dc3dinvres('step_Callback',gcbo,[],guidata(gcbo))
    case '3'
        dc3dinvres('showvol_Callback',gcbo,[],guidata(gcbo))
    case 'A'
        dc3dinvres('adddata_Callback',gcbo,[],guidata(gcbo))
    case 'C'
        dc3dinvres('showCoverage_Callback',gcbo,[],guidata(gcbo))
    case 'D'
        dc3dinvres('showdata_Callback',gcbo,[],guidata(gcbo))
    case 'E'
        dc3dinvres('exportmod_Callback',gcbo,[],guidata(gcbo))
    case 'F'
        dc3dinvres('dataforward_Callback',gcbo,[],guidata(gcbo))
    case 'G'
        dc3dinvres('optmal_Callback',gcbo,[],guidata(gcbo))
    case 'H'
        dc3dinvres('reset_Callback',gcbo,[],guidata(gcbo))
    case 'I'
        dc3dinvres('inversion_Callback',gcbo,[],guidata(gcbo))
    case 'K'
        dc3dinvres('showlcurve_Callback',gcbo,[],guidata(gcbo))
    case 'L'
        dc3dinvres('load_Callback',gcbo,[],guidata(gcbo))
    case 'M'
        dc3dinvres('showmod_Callback',gcbo,[],guidata(gcbo))
    case 'O'
        dc3dinvres('readdata_Callback',gcbo,[],guidata(gcbo))
    case 'P'
        dc3dinvres('readpro_Callback',gcbo,[],guidata(gcbo))
    case 'Model.R'
        dc3dinvres('showresponse_Callback',gcbo,[],guidata(gcbo))
    case 'S'
        dc3dinvres('save_Callback',gcbo,[],guidata(gcbo))
    case 'T'
        dc3dinvres('timelapse_Callback',gcbo,[],guidata(gcbo))
    case 'W'
        dc3dinvres('writedata_Callback',gcbo,[],guidata(gcbo))
    case 'X'
        dc3dinvres('ende_Callback',gcbo,[],guidata(gcbo))
end


% --------------------------------------------------------------------
function exportmod_Callback(hObject, eventdata, handles)
global datfile INV libmmfile
outfile=datfile;
if ishandle(1), 
    [outfile,iseps,ispdf]=getgrfile(datfile);
    if ~isempty(outfile),
        messg(['Exporting model figure to image ' outfile]);
        figure(1);
        if iseps, % Testversion %n. Zeile einklammern
            if ~isequal(libmmfile,4), errordlg('EPS export not available in Test version');return; end
            epsprint(1,outfile,ispdf);
        else exportpng(1,outfile); end
        ch=get(1,'Children');mima=get(ch(end),'Clim');%caxis;
        if INV.lolo, mima=10.^mima; end
        figure(8);cbar(mima(1),mima(2),INV.lolo,0,9,0,'\rho in \Omega m');
        if iseps, 
            epsprint(8,strrep(outfile,'.eps','-cbar.eps'),ispdf);
        else exportpng(8,strrep(outfile,'.png','-cbar.png')); end
        if ishandle(8), close(8); end
        % write text file with options and logs
        fid=fopen([outfile '.txt'],'w');
        if fid>=0,
            ss=get(handles.output,'String');
            fprintf(fid,'%s\n',ss{:});
            fprintf(fid,'Inversion Options\n');
            smethod={'Matrix','TSVD','SIRT'};
            sweight={'Equal','Smooth1','Smooth2','Coverage','resolution'};
            sauto={'Manual','L-curve','Fixed'};
            sredu={'none','delete','combine'};
            fprintf(fid,'%s inversion,weight=%s,auto=%s,reduction=%s,lambda=%.1f,keep=%s\n',...
                smethod{INV.method+1},sweight{INV.weight+1},sauto{INV.auto+1},...
                sredu{INV.redu+1},INV.lam,sauto{INV.const+1});
            %fprintf(fid,'Forward Options\n');
            fclose(fid);
        end  % text file
    else %not empty
        outfile=datfile;
    end % not empty
end
[pname,fname,ext]=fileparts(outfile);
if ishandle(2),     
    figure(2);
    [doutfile,iseps,ispdf]=getgrfile(strrep(outfile,ext,['-data' ext]));
    if ~isempty(doutfile),
        if iseps, epsprint(2,doutfile,ispdf);            
        else exportpng(2,doutfile); end
        messg(['Exporting data figure to image ' doutfile]);
    end
end
if ishandle(3),     
    figure(3);
    [routfile,iseps,ispdf]=getgrfile(strrep(outfile,ext,['-response' ext]));
    if ~isempty(routfile),
        if iseps, epsprint(3,routfile);            
        else exportpng(3,routfile); end
        messg(['Exporting response figure to image ' routfile]);
    end
end
if ishandle(5),
    figure(5);     
    [routfile,iseps,ispdf]=getgrfile(strrep(outfile,ext,['-3d' ext]));
    if ~isempty(routfile),
        if iseps, % Testversion %n. Zeile einklammern
            global libmmfile
            if ~isequal(libmmfile,4), errordlg('EPS export not available in Test version');return; end
            epsprint(5,routfile);            
        else exportpng(5,routfile); end
        messg(['Exporting 3d figure to image ' routfile]);
    end
end


% --------------------------------------------------------------------
function writedata_Callback(hObject, eventdata, handles)
global N datfile
[path,name,ext]=fileparts(datfile);
outfile=strrep(datfile,ext,'.dat');
if isempty(outfile), outfile='*.dat'; end
[fname,pname]=uiputfile(outfile,'Save Datum Points');
if fname~=0,
    datfile=strrep(fullfile(pname,fname),[pwd filesep],'');
    [ff,pp,ee]=fileparts(datfile);
    if strcmp(ee,''), datfile=[datfile '.dat']; end
    messg(['Saving data to file ',datfile]);
    saveinv3dfile(datfile,N);
    datfile=strrep(datfile,[pwd filesep],''); 
    try 
        set(hObject,'Name',[ 'DC3DInvRes - ' datfile ]);
    catch
        set(handles.figure1,'Name',[ 'DC3DInvRes - ' datfile ]);
    end
end

% --------------------------------------------------------------------
function savelb_Callback(hObject, eventdata, handles)
global N datfile
[path,name,ext]=fileparts(datfile);
outfile=strrep(datfile,ext,'.dat');
[fname,pname]=uiputfile(outfile,'Save Datum Points');
if fname~=0,
    datfile=strrep(fullfile(pname,fname),[pwd filesep],'');
    [ff,pp,ee]=fileparts(datfile);
    if strcmp(ee,''), datfile=[datfile '.dat']; end
    messg(['Exporting to Res2dinv file ',datfile]);
    saveres3dinvfile(datfile,N);
    try 
        set(hObject,'Name',[ 'DC3DInvRes - ' datfile ]);
    catch
        set(handles.figure1,'Name',[ 'DC3DInvRes - ' datfile ]);
    end
end

% --------------------------------------------------------------------
function errorest_Callback(hObject, eventdata, handles)
daterr;

% --------------------------------------------------------------------
function modres_Callback(h, eventdata, handles, varargin)
global Model N MAL s VD VM S RM RD
if isempty(RM),
    dc3dinvres('compsvd_Callback',gcbo,[],guidata(gcbo))
end
Rm=modelupdate(Model,diag(RM),0);
mal=MAL;
mal.cauto=0;mal.clog=0;mal.cmin=0;mal.cmax=1;
set(figure(1),'MenuBar','none','NumberTitle','off','Name','Model Resolution');
iconify(1);
[cmin,cmax]=draw3dmodel(Rm,MAL);
malcbar(handles,cmin,cmax,mal.clog,'res');

% --------------------------------------------------------------------
function varargout = datres_Callback(h, eventdata, handles, varargin)
global Model N MAL s S RM RD
if isempty(RD),
    dc3dinvres('compsvd_Callback',gcbo,[],guidata(gcbo))
end
mal=MAL;
mal.cauto=0;mal.clog=0;mal.cmin=0;mal.cmax=1;
Rd=diag(RD);
set(figure(2),'MenuBar','none','NumberTitle','off','Name','Data Resolution');
iconify(2);
[cmin,cmax]=plotprofiles(N,Rd,mal);
malcbar(handles,cmin,cmax,mal.clog,'res');

% --------------------------------------------------------------------
function clearsvd_Callback(h, eventdata, handles, varargin)
global s RM RD N datfile Model INV
RM=[];RD=[];s=[];
clear global RM RD s

% --------------------------------------------------------------------
function compsvd_Callback(h, eventdata, handles, varargin)
global s RM RD N datfile Model INV
[path,name,ext]=fileparts(datfile);
if iscell(Model.M),
    matfile=strrep(datfile,ext,'-s.mat');
else
    matfile=strrep(datfile,ext,'-sens.mat');
end
if exist(matfile,'file'),
    load(matfile);
    if (~isequal(Model.z,zsave)),
        global S
        dc3dinvres('calcsens_Callback',gcbo,[],guidata(gcbo));
    end
else
    global S 
    dc3dinvres('calcsens_Callback',gcbo,[],guidata(gcbo))
end
t0=clock;
data=length(N.r);
if 1,
    if ~isfield(N,'err'),
        daterr;
        messg(sprintf('Error min=%.1f%% max=%.1f%% mean=%.1f%%',...
            min(err)*100,max(err)*100,mean(err)*100));
    end
    D=spdiags(1./log(1+N.err),0,data,data);
    messg('Computing singular value decomposition...');
    set(gcf,'Pointer','watch');
    [VD,s,VM]=svd(D*S);  % Model & Data Vectors
    set(gcf,'Pointer','arrow');
    s=diag(s);
    %dR=log(1+randn(size(N.r)).*err);
    dR=log(N.r)-log(Model.Bg(1));
    lambda=0;
    if isfield(INV,'first')&&(INV.first>0),
        lambda=INV.first;
    elseif isfield(INV,'lam')&&(INV.lam>0),
        lambda=INV.lam;
    else
        %[dM,lambda]=invshifted(D*S,D*dR,100,0.001,0.5);
        lambda=10;
    end
else % not used anymore
    [VD,s,VM]=svd(S);  % Model & Data Vectors
    s=diag(s);
    lambda=0.05;
end
messg(sprintf('ready(%.2fs) min(sv)= %g, max(sv)= %g',...
    etime(clock,t0),min(s),max(s)));
% Regularization
r=max(find(s));s=s(1:r);f=s.^2./(s.^2+lambda);
RD=VD(:,1:r)*diag(f)*VD(:,1:r)';
Rd=diag(RD);
RM=VM(:,1:r)*diag(f)*VM(:,1:r)';
messg(sprintf('Inf. content %.1f, effectivity %d%%, model inf. %d%%',...
    sum(Rd),round(mean(Rd)*100),round(mean(diag(RM))*100)));


% --------------------------------------------------------------------
function setrasd_Callback(hObject, eventdata, handles)
global N Model
N.r=Model.R;
if ishandle(2),
    dc3dinvres('showdata_Callback',gcbo,[],guidata(gcbo))
end

% --------------------------------------------------------------------
function noisifyr_Callback(hObject, eventdata, handles)
global N Model
daterr;
noise=randn(size(N.r)).*N.err;
Model.R=Model.R.*(1+noise);
Model.isfor=1;

% --- Executes during object creation, after setting all properties.
function slidesens_CreateFcn(hObject, eventdata, handles)
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function slidesens_Callback(hObject, eventdata, handles)
global S MAL Model N INV malstat
D=spdiags(1./log(1+N.err(:)),0,length(N.err),length(N.err));
nr=round(get(handles.slidesens,'Value'));
if (~isempty(malstat))&&(malstat>0), % cell resolution
    MCR=modelcellres(S,INV.lam,getcmatrix3d(Model,INV),D,nr);
    if malstat==3,
        set(figure(5),'MenuBar','none','NumberTitle','off','Name','3D-Resolution');
        iconify(5); % früher m v o l
        modelvolume(reshape(MCR,size(Model.M)),Model.x,Model.y,Model.z,max(MCR(:))/3,'.','.','.');
    else
        mm=max(abs(MCR(:)));
        mal=struct('cmap',2,'cauto',0,'cmin',-mm,'cmax',mm);
        set(figure(1),'MenuBar','none','NumberTitle','off','Name','Cell Resolution');
        iconify(1);
        draw3dmodel(modelupdate(Model,MCR,0),mal)        
    end
else
    set(figure(1),'MenuBar','none','NumberTitle','off','Name','Sensitivity');
    iconify(1);
    mal=MAL;mal.clog=0;mal.cmap=2;mal.cauto=0;mal.elec=0;
    mal.cmax=0.01;mal.cmin=-mal.cmax;mal.alpha=0;
    if nr<1, nr=1; end
    if nr>size(S,1), nr=size(S,1); end
    ss=sprintf('%d: A(%g,%g)',nr,N.elec(N.a(nr),1),N.elec(N.a(nr),2));
    if N.b(nr)>0,
        ss=sprintf('%s B(%g,%g)',ss,N.elec(N.b(nr),1),N.elec(N.b(nr),2));
    end
    ss=sprintf('%s M(%g,%g)',ss,N.elec(N.m(nr),1),N.elec(N.m(nr),2));
    if N.n(nr)>0,
        ss=sprintf('%s N(%g,%g)',ss,N.elec(N.n(nr),1),N.elec(N.n(nr),2));
    end
    set(handles.status,'String',ss);
    draw3dmodel(modelupdate(Model,full(S(nr,:)),0),mal);
    malcbar(handles,mal.cmin,mal.cmax,mal.clog,'');

%     elecs=N.elec([N.a(nr) N.b(nr) N.m(nr) N.n(nr)],1:2)
%     ch=findobj(get(1,'Children'),'Type','Axes');axes(ch(end));
%     hold on;plot(elecs(:,1),elecs(:,2),'k.');hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = showsens_Callback(h,eventdata,handles)
global S malstat
malstat=0;
set(handles.slidesens,'Visible','On','Value',0);
ma=size(S,1);
set(handles.slidesens,'sliderstep',[1 1]*(1/(ma-1)),'max',ma,'min',1,'Value',1);
dc3dinvres('slidesens_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function rmsout_Callback(hObject, eventdata, handles)
global RMS CHIQ
% RMS=[RMS rms(N.r,Model.R,INV.lolo)];
% CHIQ=[CHIQ chi2(N.r,Model.R,INV.lolo)];
messg(['CHI^2 = ',sprintf('%.2f ',CHIQ) '(RMS=' sprintf('%.2f%%',RMS(end)) ')']);


% --------------------------------------------------------------------
function display_Callback(hObject, eventdata, handles)
global Model N S Coverage
chon='Off';ipon='Off';
if isfield(Model,'IP')&&(isequal(size(Model.M),size(Model.IP))), chon='On'; end
if isfield(N,'ip')&&(isequal(size(N.ip),size(N.r))), ipon='On'; end
set(handles.showip,'Enable',ipon);
set(handles.showcharge,'Enable',chon);
set(handles.showmod,'Enable',ooemp(Model));
set(handles.showvol,'Enable',ooemp(Model));
set(handles.showsens,'Enable',ooemp(Model,S));
set(handles.showcov,'Enable',ooemp(Model,Coverage));
oodata=ooemp(N);
set(handles.showdata,'Enable',oodata);
set(handles.showresponse,'Enable',oodata);
set(handles.showmisfit,'Enable',oodata);
set(handles.showerror,'Enable',oodata);
    
% --------------------------------------------------------------------
function showip_Callback(hObject, eventdata, handles)
global N MAL
set(figure(2),'MenuBar','none','NumberTitle','off','Name','Apparent Chargeability');
iconify(2);
mal=MAL;mal.vie=1;mal.clog=0;mal.cauto=1;mal.canot='\phi in mrad';
[cmin,cmax]=plotprofiles(N,N.ip,mal);
malcbar(handles,cmin,cmax,0,'m_a[mrad]');
set(handles.status,'String','Ready');
set(handles.slidesens,'Visible','Off');

% --------------------------------------------------------------------
function showcharge_Callback(hObject, eventdata, handles)
global MAL S Model
set(figure(1),'MenuBar','none','NumberTitle','off','Name','Chargeability');
iconify(1);
mal=MAL;mal.clog=0;mal.cauto=1;mal.alpha=0;
% mal.cauto=0;mm=max(abs(Model.IP(:)));mal.cmin=-mm;mal.cmax=mm;mal.cmap=2;
[cmin,cmax]=draw3dmodel(modelupdate(Model,Model.IP,0),mal);
malcbar(handles,cmin,cmax,0,'mrad');
set(handles.status,'String','Ready');
set(handles.slidesens,'Visible','Off');

% --------------------------------------------------------------------
function ipinv_Callback(hObject, eventdata, handles)
global Model N INV S
inv=INV;inv.auto=2;
Model.IP=zeros(size(Model.M));
fi=find(N.ip>0);
Model.IP(:)=minvers(S(fi,:),N.ip(fi),inv,N.err); 
% Model.IP(:)=minvers(S,N.ip,inv,N.err,Model);
dc3dinvres('showcharge_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function data_Callback(hObject, eventdata, handles)
global N Model Model
oodata=ooemp(N);
set(handles.visitdata,'Enable',oodata);
set(handles.dataedit,'Enable',oodata);
set(handles.datashowel,'Enable',oodata);
set(handles.datadeldead,'Enable',oodata);
set(handles.writedata,'Enable',oodata);
set(handles.savelb,'Enable',oodata);
set(handles.datasavepro,'Enable',oodata);
set(handles.dataforward,'Enable',ooemp(Model,N));
if isfield(Model,'R'),
    set(handles.setrasd,'Enable',ooemp(Model.R,N));
    set(handles.noisifyr,'Enable',ooemp(Model.R,N));
else
    set(handles.setrasd,'Enable','Off');
    set(handles.noisifyr,'Enable','Off');
end

% --------------------------------------------------------------------
function dataedit_Callback(hObject, eventdata, handles)
global N Model INV CHIQ RMS
ndelete;
RMS=rms(N.r,Model.R,INV.lolo);
CHIQ=chi2(N.r,Model.R,N.err,INV.lolo);
dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function datasavepro_Callback(hObject, eventdata, handles)
global N datfile
[path,name,ext]=fileparts(datfile);
outfile=strrep(datfile,ext,'.pro');
if isempty(outfile), outfile='*.pro'; end
[fname,pname]=uiputfile(outfile,'Save Datum Points');
if fname~=0,
    profile=strrep(fullfile(pname,fname),[pwd filesep],'');
    [ff,pp,ee]=fileparts(profile);
    if strcmp(ee,''), profile=[profile '.pro']; end
    messg(['Saving profile data to file ',profile]);
    savepro(profile,N);
end

% --------------------------------------------------------------------
function datashowel_Callback(hObject, eventdata, handles)
global N
if ~isfield(N,'elec')||(ndims(N.elec)<2), return; end
set(figure(2),'MenuBar','none','NumberTitle','off','Name','Electrode Positions');
iconify(2);
clf;
plot(N.elec(:,1),N.elec(:,2),'+');
grid on


% --------------------------------------------------------------------
function datadeldead_Callback(hObject, eventdata, handles)
global N
N=deldeadelecs(N);
if isfield(N,'zweid'),
    for i=1:length(N.zweid), N.zweid{i}=deldeadelecs(N.zweid{i}); end
end
dc3dinvres('showdata_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function dataforward_Callback(hObject, eventdata, handles)
global N Model RMS CHIQ INV FOR
Model.R=mfdfwd3d(Model,N,FOR);
Model.isfor=1;
RMS=[RMS rms(N.r,Model.R,INV.lolo)];
CHIQ=[CHIQ chi2(N.r,Model.R,N.err,INV.lolo)];
dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))


% --------------------------------------------------------------------
function exportcbar_Callback(hObject, eventdata, handles)
global datfile
[path,name,ext]=fileparts(datfile); 
if isempty(ext), outfile='cbar.png';
else outfile=strrep(datfile,ext,'-vcbar.png'); end
[fname,pname]=uiputfile({'*.png ';'*.eps ';'*.pdf'},'Save Figure as',outfile);
iseps=strfind(fname,'.eps');
ispdf=strfind(fname,'.pdf');
if fname~=0,
    [pp,nn,ee]=fileparts(fname);
    nn=strrep(nn,'.','_');
    pname=strrep(pname,[pwd filesep],'');
    outfile=fullfile(pname,[nn ee]);
    if ispdf, iseps=1;outfile=strrep(outfile,'.pdf','.eps'); end
    set(figure(8),'MenuBar','none');
    copyobj(handles.cbar,8);
    if iseps, 
        epsprint(8,strrep(outfile,'.eps',''),~isempty(ispdf));
        if ispdf, delete(outfile); end
    else exportpng(8,outfile); end
    close(8);
end


% --------------------------------------------------------------------
function cellres_Callback(hObject, eventdata, handles)
global S malstat
malstat=2;
ma=size(S,2);
set(handles.slidesens,'Visible','On','Value',fix(ma/2)-5);
set(handles.slidesens,'sliderstep',[1 1]*(1/(ma-1)),'max',ma,'min',1,'Value',fix(ma/2));
dc3dinvres('slidesens_Callback',gcbo,[],guidata(gcbo))

% --------------------------------------------------------------------
function cellres3d_Callback(hObject, eventdata, handles)
global S malstat
malstat=3;
ma=size(S,2);
set(handles.slidesens,'Visible','On','Value',fix(ma/2)-5);
set(handles.slidesens,'sliderstep',[1 1]*(1/(ma-1)),'max',ma,'min',1,'Value',fix(ma/2));
dc3dinvres('slidesens_Callback',gcbo,[],guidata(gcbo))
view(3);


% --------------------------------------------------------------------
function helpmenu_Callback(hObject, eventdata, handles)
% hObject    handle to helpmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function docu_Callback(hObject, eventdata, handles)
ibrowse(['file://' pwd filesep 'doc' filesep 'index.html']); 

% --------------------------------------------------------------------
function onlinedoc_Callback(hObject, eventdata, handles)
ibrowse(['http://dc3dinvres.resistivity.net']);

% --------------------------------------------------------------------
function webinv_Callback(hObject, eventdata, handles)
% ibrowse(['http://red-bull.geophysik.tu-freiberg.de/~guenther/webinv/index3d.html']);
ibrowse(['http://resistivity.net/index.php?webinv']); 

% --------------------------------------------------------------------
function exportvtk_Callback(hObject, eventdata, handles)
global Model datfile Coverage
[path,name,ext]=fileparts(datfile);
outfile=strrep(datfile,ext,'.vtk');
[fname,pname]=uiputfile(outfile,'Export model to VTK-File');
if fname==0, return; end% no valid name
outfile=fullfile(pname,fname);
[ff,pp,ee]=fileparts(outfile);
if strcmp(ee,''), outfile=[outfile '.vtk']; end
messg(['Exporting model to VTK File : ' outfile]);
% errordlg('VTK export not available in Test version');return;
vtkexport3d(outfile,Model);

% --------------------------------------------------------------------
function reweighting_Callback(hObject, eventdata, handles)
global N Model CHIQ INV RMS CHIQ
dR=(log(N.r)-log(Model.R))./log(1+N.err);
w=abs(dR)*sum(abs(dR))/sum(dR.^2);
w(find(w<1))=1;
N.err=N.err.*w;
RMS=[RMS rms(N.r,Model.R,INV.lolo)];
CHIQ=[CHIQ chi2(N.r,Model.R,N.err,INV.lolo)];
dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
set(figure(2),'MenuBar','none','NumberTitle','off','Name','Reweighting');
plotprofiles(N,w);

% --------------------------------------------------------------------
function visitdata_Callback(hObject, eventdata, handles)
viewprofiles;

% --------------------------------------------------------------------
function filemenu_Callback(hObject, eventdata, handles)
global N
set(handles.adddata,'Enable',ooemp(N));

% --------------------------------------------------------------------
function modelmenu_Callback(hObject, eventdata, handles)
global Model
oomod=ooemp(Model);
set(handles.optmod,'Enable',oomod);
set(handles.modfree,'Enable',oomod);
set(handles.modexp,'Enable',oomod);
set(handles.modimp,'Enable',oomod);
set(handles.exportvtk,'Enable',oomod);

% --------------------------------------------------------------------
function options_Callback(hObject, eventdata, handles)
global N
oodata=ooemp(N);
set(handles.errorest,'Enable',oodata);
set(handles.reweighting,'Enable',oodata);

% --------------------------------------------------------------------
function invmenu_Callback(hObject, eventdata, handles)
global Model N S
oomod=ooemp(Model);
oomoddat=ooemp(Model,N);
set(handles.reset,'Enable',oomod);
set(handles.start,'Enable',oomoddat);
set(handles.step,'Enable',oomoddat);
set(handles.inversion,'Enable',oomoddat);
set(handles.calcsens,'Enable',oomoddat);
set(handles.delsens,'Enable',ooemp(S));
set(handles.ipinv,'Enable',onoff(~isempty(N)&&isfield(N,'ip')&&~isempty(S)));

% --------------------------------------------------------------------
function resmenu_Callback(hObject, eventdata, handles)
global Model S RM
oores=ooemp(S);
set(handles.cellres,'Enable',oores);
set(handles.cellres3d,'Enable',oores);
set(handles.compsvd,'Enable',oores);
oorm=ooemp(RM);
set(handles.modres,'Enable',oorm);
set(handles.datres,'Enable',oorm);
set(handles.clearsvd,'Enable',oorm);

function res=ooemp(varargin)
res='On';
for i=1:nargin,
    if isempty(varargin{i}), res='Off'; end
end

function res=onoff(logi)
res='Off';
if ~isempty(logi)&&(logi), res='On'; end


% --------------------------------------------------------------------
function setasreference_Callback(hObject, eventdata, handles)
global Model Mref
if iscell(Model),
  Mref=[];for k=1:length(Model.M), Mref=[Mref;Model.M{k}(:)]; end
else Mref=Model.M(:); end

% --------------------------------------------------------------------
function exportohm_Callback(hObject, eventdata, handles)
global datfile N
[pname,fname,ext]=fileparts(datfile);
outfile=strrep(datfile,ext,'.ohm');
[fname,pname]=uiputfile(outfile,'Save Ohm-File');
if isequal(fname,0), return; end
ohmfile=strrep(fullfile(pname,fname),[pwd filesep],'');
[pname,fname,ee]=fileparts(ohmfile);
if strcmp(ee,''), ohmfile=[ohmfile '.ohm']; end
messg(['Exporting to ohm-xz file ',ohmfile]);
if isfield(N,'topo'),
    T=rmfield(N,'topo');
    if size(N.elec,1)==size(N.topo), T.elec=N.topo; end %simple case/all have topo
%     T.elec(:,3)=interp2(N.topo(:,1),
else
    T=N;
end
T.rho=N.r./N.k;
if isfield(T,'r'), T=rmfield(T,'r'); end
if isfield(T,'i'), T=rmfield(T,'i'); end
if isfield(T,'u'), T=rmfield(T,'u'); end
if isfield(T,'ip'), T=rmfield(T,'ip'); end
saveinv2dfile(ohmfile,T);

% --- Executes on button press in full.
function full_Callback(hObject, eventdata, handles)
dc3dinvres('inversion_Callback',gcbo,[],guidata(gcbo))

% --- Executes on button press in full.
function timelapse_Callback(hObject, eventdata, handles)
global x z Model INV datfile MAL S N
if iscell(Model.M), return; end
[pname,fname,ext]=fileparts(datfile);
infile=strrep(datfile,[fname ext],'*.dat');
[newfile,newpath]=uigetfile(infile,'Time lapse file!');
if isequal(newfile,0), return; end
N1=read3dfile(fullfile(newpath,newfile));
dM=Model.M;
if ~isequal(N.elec,N1.elec), uiwait(errordlg('Electrodes not identical!','Time lapse error'));return; end
% [c,ia,ib]=intersect([N1.a N1.b N1.m N1.n],[N.a N.b N.m N.n],'rows');
[ggN,ia,ib]=dcintersect(N1,N);
if ~any(ia), messg('no intersecting data!');return; end
saveinv3dfile(fullfile(newpath,'ratio.dat'),ggN);
messg(sprintf('Found %d intersecting data and wrote ratio.dat',length(ib)));
if INV.lolo, dR=log(N1.r(ia))-log(N.r(ib));
else dR=N1.r(ia)-N.r(ib); end
err=N.err(ib);
if isfield(N1,'err'), err=(N1.err(ia)+N.err(ib))/2; end
dM(:)=exp(minvers(S(ib,:),dR,INV,err,Model));
% ggN.elec=N.elec;ggN.a=N.a(ib);ggN.b=N.b(ib);ggN.m=N.m(ib);
% ggN.n=N.n(ib);ggN.err=(N1.err(ia)+N.err(ib))/2;ggN.r=N1.r(ia)./N.r(ib);
mal=MAL;mal.cauto=0;
if INV.lolo,
    mm=max(max(dM(:)),1/min(dM(:)));
    mal.clog=1;mal.cmin=1./mm;mal.cmax=mm;
else
    mm=max(abs(dM));
    mal.clog=0;mal.cmin=-mm;mal.cmax=mm;
end
figure(1);
[cmin,cmax]=draw3dmodel(modelupdate(Model,dM,0),mal);
malcbar(handles,cmin,cmax,mal.clog,'\Delta in %');


% --------------------------------------------------------------------
function electrodevtk_Callback(hObject, eventdata, handles)
global N datfile
[path,name,ext]=fileparts(datfile);
outfile=strrep(datfile,ext,'-elec.vtk');
[fname,pname]=uiputfile(outfile,'Export electrodes to VTK-File');
if fname==0, return; end
outfile=fullfile(pname,fname);
[ff,pp,ee]=fileparts(outfile);
if strcmp(ee,''), outfile=[outfile '.vtk']; end
savepointvtk(outfile,N.elec);


% --------------------------------------------------------------------
function comparemodel_Callback(hObject, eventdata, handles)
global Model datfile MAL
[fpath,name,ext]=fileparts(datfile);
if isempty(ext),
    outfile='*.mod';
else
    outfile=strrep(datfile,ext,'.mod');
end
[fname,pname]=uigetfile(outfile,'Model to compare');
if fname==0, return; end

[Mod1.M,Mod1.x,Mod1.y,Mod1.z,Coverage]=modelimport3d(fullfile(pname,fname));
for i=1:size(Mod1.M,3), mm=Mod1.M(:,:,1);Mod1.Bg(i)=median(mm(:)); end
Mod1.Bg(end+1)=Mod1.Bg(end);
if isequal(Model.x(:),Mod1.x(:))&&isequal(Model.y(:),Mod1.y)&&isequal(Model.z(:),Mod1.z(:)),
    Mod1.M=Model.M./Mod1.M;
    mm=max(max(Mod1.M(:)),1/min(Mod1.M(:)));
    mal=MAL;mal.cauto=0;mal.clog=1;mal.cmin=1/mm;mal.cmax=mm;mal.cmap=2;
    figure(1);draw3dmodel(Mod1,mal);
    malcbar(handles,mal.cmin,mal.cmax,mal.clog,'');
    %%
    ss=sprintf('Model Difference RMS = %.2f%%',rms(Model.M(:),Mod1.M(:))); 
    if strcmp(questdlg('Export ratio as model file?',ss),'Yes'),
        [fpath,name,ext]=fileparts(datfile);
        [fname,pname]=uiputfile(fullfile(fpath,'ratio.mod'),'Export model to ASCII-File');
        if fname~=0,
            outfile=fullfile(pname,fname);
            [fpath,name,ee]=fileparts(outfile);
            if strcmp(ee,''), outfile=[outfile '.mod']; end
%             sec=Mod.Cov;isip=0;
%             if isfield(Model,'IP')&&isequal(size(Mod.M),size(Mod.IP)),
%                 sec=Mod.IP;isip=1; end
            modelexport3d(outfile,Mod1,Model.x,Model.y,Model.z,Coverage);
        end
    end % save as ratio.mod
else
    errordlg('Model parameterization is not identical!');    
end


% --------------------------------------------------------------------
function about_Callback(hObject, eventdata, handles)
f=msgbox({'DC3dInvRes';'Version 2.11.6';'Author: Thomas Günther';'Resistivity.net'},'About');
iconify(f);uiwait(f);
