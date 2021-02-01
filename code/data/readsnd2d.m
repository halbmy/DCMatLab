function N=readsnd2d(filename)

% READSND2D - Read 2d file by use of soundings
% N = readsnd2d(filename);
% Format:
% sounding_filename_1 position
% ...
% sounding files consists of 3 columns for AB/2 MN/2 and apparent resistivity

N=[];
[fpath,fname,fext]=fileparts(filename);
fid=fopen(fullfile(fpath,[fname fext]),'r');
if isequal(fid,-1), error('Could not open file!'); end
try
    zeile='';
    while isempty(zeile), zeile=destrip(fgetl(fid)); end
    l=0;
    while ischar(zeile),
        l=l+1;
        vesfile=sscanf(zeile,'%s %*s %*s');
        pos=sscanf(zeile,'%*s %f');
        type=sscanf(zeile,'%*s %*s %d');        
        fprintf('Sounding %s at x=%.1fm\n',vesfile,pos);
        vesfull=fullfile(fpath,vesfile);
        if ~exist(vesfull,'file'),
            fprintf('VES file %s does not exist!\n',vesfile);
        else
            data=readvesfile(vesfull);
            AA=-data(:,1);MM=-data(:,2);NN=-MM;BB=-AA;
            V=abmn2n(AA,BB,MM,NN);
            if type==2, V.b(:)=0; end
            if type==1, V.a=V.b;V.b(:)=0;du=V.m;V.m=V.n;V.n=du; end
            V.elec(:,1)=V.elec(:,1)+pos;
            V.elec(:,2)=0;
            V.r=data(:,3);
            numbers{l}=(1:length(V.r));
            if l>1, numbers{l}=numbers{l}+length(N.r); end
            if l==1, N=V; else N=combdata2d(N,V); end
            zeile='';
            while isempty(zeile), zeile=destrip(fgetl(fid)); end
            eind{l}=[data(:,1) data(:,3)];
            names{l}=strrep(vesfile,'.ves','');
            positions(l,1)=pos;
        end
    end
catch
    display(lasterr);
end %try
fclose(fid);
if ~isempty(N), N.k=getkonf(N);N=deldeadelecs(N); end
if exist('eind','var'), N.eind=eind; end
if exist('positions','var'), N.pos=positions; end
if exist('names','var'), N.names=names; end
if exist('numbers','var'), N.nr=numbers; end
% 
% function N=abmn2n(AA,BB,MM,NN)
% [N.elec,SI,SJ]=unique([AA;BB;MM;NN],'rows');
% [TF,LOC]=ismember(AA,N.elec);
% N.a=LOC;
% [TF,LOC]=ismember(BB,N.elec);
% N.b=LOC;
% [TF,LOC]=ismember(MM,N.elec);
% N.m=LOC;
% [TF,LOC]=ismember(NN,N.elec);
% N.n=LOC;
% 
% function N=combdata2d(N,N1)
% 
% data=length(N.a);
% ne=size(N.elec,1);
% ne1=size(N.elec,1);
% data1=length(N1.a);
% 
% index=(1:size(N1.elec,1))'+ne;
% %Elektroden anhängen
% [aa,bb]=meshgrid(N1.elec(:,1)+N1.elec(:,2)*12.34,N.elec(:,1)+N.elec(:,2)*12.34);
% [ii,jj]=find((aa-bb)==0);
% index(jj)=ii;
% ind=find(index>ne);
% N.elec=[N.elec;N1.elec(ind,:)];
% cu=cumsum(index>ne);
% index(ind)=ne+cu(ind);
% N.a(data+1:data+data1)=index(N1.a);
% N.b(data+data1)=0; % verlängern
% fb=find(N1.b>0);
% N.b(fb+data)=index(N1.b(fb));
% N.m(data+1:data+data1)=index(N1.m);
% N.n(data+data1)=0; % verlängern
% fn=find(N1.n>0);
% N.n(fn+data)=index(N1.n(fn));
% if isfield(N,'r')&&isfield(N1,'r'), N.r=[N.r(:);N1.r(:)]; end
% if isfield(N,'k')&&isfield(N1,'k'), N.k=[N.k(:);N1.k(:)]; end
% if isfield(N,'err')&isfield(N1,'err'),
%     N.err=[N.err(:);N1.err(:)];
% else N.err=[]; end
% if isfield(N,'ip')&isfield(N1,'ip'),
%     N.ip=[N.ip(:);N1.ip(:)];
% else N.ip=[]; end
% if isfield(N,'i')&isfield(N1,'i'),
%     N.i=[N.i(:);N1.i(:)];
% else N.i=[]; end
% if isfield(N,'u')&isfield(N1,'u'),
%     N.u=[N.u(:);N1.u(:)];
% else N.u=[]; end
