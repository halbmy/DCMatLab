function [Mesh,N]=postmodel2d(filename)

% POSTMODEL2D - postprocessor for 2d BERT inversion
% postmodel(filename) 
% filename is either a cfg file (e.g. in a result directory) or a zip file
% sent back from web inversion

[fpath,fname,fext]=fileparts(filename);
if isequal(fext,'.zip'),
    tdir=tempname;
    fi=strfind(tdir,filesep);    
    mkdir(tdir(1:fi(end)-1),tdir(fi(end)+1:end));
%     mkdir(tdir);    
    unzip(filename,tdir);
    aa=dir(tdir);
    if length(aa)<4, % result dir in zip file
      tdir=[tdir filesep aa(end).name];
      aa=dir(tdir);
    end
    for i=1:length(aa), if strfind(aa(i).name,'.cfg');break; end; end
    filename=[tdir filesep aa(i).name];
end
fpath=fileparts(filename);if isempty(fpath), fpath='.'; end
meshname=[fpath filesep 'mesh' filesep 'meshParaDomain'];
fid=fopen(filename,'r');zeile=fgetl(fid);istopo=0;
while ischar(zeile),
   fi=strfind(zeile,'=');first=zeile(1:fi-1);last=zeile(fi+1:end);zeile=fgetl(fid); 
   if strcmp(first,'DATAFILE'), datfile=last; end
   if strcmp(first,'TOPOGRAPHY'), istopo=(str2num(last)>0); end
end
fclose(fid);
[dpath,dname,dext]=fileparts(datfile);
datafile=[fpath filesep strrep(datfile,dext,'.data')];
N=readinv2dfile(datafile);
if isfield(N,'topo'),
    di=rndig(sqrt(sum(diff(N.topo).^2,2)));
    if length(unique(di))>5, N.elec(:,1)=(0:size(N.elec,1)-1)'*round(median(di)); end
    N.elec(:,2)=0;
end
Mesh=loadmesh(meshname);
isbert2=exist(fullfile(fpath,'bert.log'),'file');
if isbert2,
    aa=dir([fpath filesep 'model_1*.vector']);
    if length(aa)<2, aa=dir([fpath filesep 'model_*.vector']); end
    fp=strfind(aa(end).name,'.');fs=strfind(aa(end).name,'_');
    Mesh.iter=str2num(aa(end).name(fs(1)+1:fp(1)-1));
else
    aa=dir([fpath filesep 'model_iter.1*.vector']);
    if length(aa)<2, aa=dir([fpath filesep 'model_iter.*.vector']); end
    fp=strfind(aa(end).name,'.');
    Mesh.iter=str2num(aa(end).name(fp(1)+1:fp(2)-1));
end
if isempty(aa),
    display('No model file found!');return;
else
    modname=aa(end).name;
    fid=fopen([fpath filesep modname],'r');model=mytextscan(fid,'%f');fclose(fid);
    Mesh.model=model{1};
end
covfile=fullfile(fpath,'sensCov.vector');
if isbert2, covfile=fullfile(fpath,'coverage.vector'); end
if exist(covfile,'file'),
    fidcov=fopen(covfile,'r');
    senscov=fscanf(fidcov,'%f',Inf);
    fclose(fidcov);
    if ~isbert2, senscov=log(senscov); end
    [nn,hh]=hist(senscov,50);
    nnn=cumsum(nn)/length(senscov);
    mi=hh(min(find(nnn>0.02)));
    ma=hh(max(find(nnn<0.5)));
    alfa=(senscov-mi)/(ma-mi);
    alfa(alfa<0)=0;
    alfa(alfa>1)=1;
else
    alfa=ones(size(Mesh.model));
end
Mesh.alfa=alfa;
if isempty(aa),
    display('No model response file found!');return;
else
    if isbert2,
        morname='response.vector';
    else
        aa=dir([fpath filesep 'modelResponse.1*.vector']);
        if length(aa)<2, aa=dir([fpath filesep 'modelResponse.*.vector']); end
        morname=aa(end).name;
    end
    morfile=fullfile(fpath,morname);
    if exist(morfile,'file'),
        fid=fopen(morfile,'r');response=mytextscan(fid,'%f');fclose(fid);
        response=response{1};misfit=(response./N.r-1)*100;
        N.response=response;
    end
end
ipfile=fullfile(fpath,'ip_model.vector');
if isbert2, ipfile='phase.vector'; end
if exist(ipfile,'file'),
   fid=fopen(ipfile,'r');erg=mytextscan(fid,'%f');fclose(fid);
   Mesh.ipmodel=erg{1};
end
if isequal(fext,'.zip'),
    rmdir(tdir,'s');
end
di=max(Mesh.node)-min(Mesh.node); 
figure(1);clf;set(1,'Renderer','zbuffer');
set(gcf,'MenuBar','none','NumberTitle','off','Name','Inversion result');
po=get(gcf,'Position');po(1)=20;po(2)=50;po(3)=800;po(4)=di(2)/di(1)*1.3*800+50;
while po(4)>800, po(3:4)=po(3:4)/2; end
set(gcf,'Position',po); 
tripatchmod(Mesh,Mesh.model,alfa,struct('canot','Ohmm'));
t=text(max(Mesh.node(:,1)),min(Mesh.node(:,2)),'BERT@resistivity.net');
set(t,'HorizontalAlignment','right','VerticalAlignment','bottom');
if isfield(N,'response'),
    figure(2);clf;mal=struct('cauto',0,'clog',1,'cmin',min(N.r),'cmax',max(N.r));
    set(gcf,'MenuBar','none','NumberTitle','off','Name','Data and Misfit');
    subplot(3,1,1);showdata2d(N,N.r,mal);
    xl=get(gca,'XLim');xl=xl(1)+0.05*diff(xl);
    yl=get(gca,'YLim');yl=yl(1)+0.9*diff(yl);
    text(xl,yl,'Measured Data in \Omega m');
    subplot(3,1,2);showdata2d(N,response,mal);
    text(xl,yl,'Calculated Data in \Omega m');
    subplot(3,1,3);showdata2d(N,misfit);
    text(xl,yl,'Data Misfit in %');
end
if (istopo)&isfield(N,'t'),
    if max(N.t)==0, 
        NN=N;NN.elec=[0;cumsum(sqrt(sum(diff(N.elec).^2,2)))];
        NN.elec(:,2)=0;
        NN.k=getkonf2d(NN);N.t=N.k./NN.k;
    end
    figure(3);clf;rawdata=N.r.*N.t;
    set(gcf,'MenuBar','none','NumberTitle','off','Name','Topography effect');
    subplot(3,1,1);showdata2d(N,rawdata,mal);
    xl=get(gca,'XLim');xl=xl(1)+0.05*diff(xl);
    yl=get(gca,'YLim');yl=yl(1)+0.9*diff(yl);
    text(xl,yl,'Raw Data in \Omega m');
    subplot(3,1,2);showdata2d(N,N.t);
    text(xl,yl,'Topography effect');
    subplot(3,1,3);showdata2d(N,N.r,mal);
    text(xl,yl,'Corrected Data in \Omega m');
end
if isfield(Mesh,'ipmodel')&&(length(Mesh.ipmodel)==length(Mesh.model)),
    figure(4);clf;set(4,'Renderer','OpenGL');
    set(gcf,'MenuBar','none','NumberTitle','off','Name','IP model');
    po=get(gcf,'Position');po(1)=50;po(3)=800;po(4)=di(2)/di(1)*1.3*800+50;
    while po(4)>800, po(3:4)=po(3:4)/2; end
    set(gcf,'Position',po);
    tripatchmod(Mesh,Mesh.ipmodel,alfa);
    t=text(max(Mesh.node(:,1)),min(Mesh.node(:,2)),'BERT@resistivity.net');
    set(t,'HorizontalAlignment','right','VerticalAlignment','bottom');
end
