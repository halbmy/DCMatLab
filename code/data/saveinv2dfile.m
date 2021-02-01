function saveinv2dfile(fname,N,maxi,eltoken)

% SAVEINV2DFILE - Save INV2D file format
% saveinv2dfile(filename,N)
% saveinv2dfile(filename,N,verbose_format)
% saveinv2dfile(filename,N,verbose_format,electrode_token)
% help readinv2dfile for format

if nargin<2, error('saveinv2dfile(filename,Data)'); end
if nargin<3, maxi=1; end
if nargin<4, eltoken='# x z'; end
newline='\n';
if maxi, newline='\r\n'; end
fid=fopen(fname,'w');
if fid<0,
    error('File not found!');
end
if maxi==0, %sorting of electrodes
  N2=N;
  [N.elec,ind]=sortrows(N2.elec,1);
  N.a=ind(N2.a);N.b=zeros(size(N.a));
  fi=find(N2.b);N.b(fi)=ind(N2.b(fi));
  N.m=ind(N2.m);N.n=zeros(size(N.a));
  fi=find(N2.n);N.n(fi)=ind(N2.n(fi));
end
fprintf(fid,'%d',size(N.elec,1));
if maxi, fprintf(fid,'# Number of electrodes'); end
fprintf(fid,newline);
if isfield(N,'topo'),
%     N.elec(:,3)=interp1(N.topo(:,1),N.topo(:,2),N.elec(:,1),'linear','extrap');
    if maxi==0, N.elec(:,2)=[]; end
    if maxi,
      fprintf(fid,[eltoken newline]); end
else
    if maxi, fprintf(fid,[eltoken newline]); end
end
ss='';
for l=1:size(N.elec,2), ss=[ss '\t%g']; end
ss=[ss newline];
ss(1:2)=[];
fprintf(fid,ss,N.elec');
fprintf(fid,'%d',length(N.a));
if maxi, fprintf(fid,'# Number of data'); end
fprintf(fid,newline);
% if maxi, fprintf(fid,['# Electrode numbers A B M N (0=inf), rhoa, error, ...' newline]); end
mess=double([N.a(:) N.b(:) N.m(:) N.n(:)]);
elst='%4d';sep='\t';%sep=' ';
ss=[elst sep elst sep elst sep elst];%sep='  ';
dtab='\t';doku=['#a' dtab 'b' dtab 'm' dtab 'n'];
if isfield(N,'r')&&(length(N.r)==length(N.a))&&maxi,
    mess=[mess N.r(:)];
    doku=[doku dtab 'rhoa'];
    if min(N.r)>1, ss=[ss sep '%.2f']; else ss=[ss sep '%g']; end
end
if isfield(N,'rho')&&(length(N.rho)==length(N.a))&&maxi,
    mess=[mess N.rho(:)];
    if min(abs(N.rho))>1, ss=[ss sep '%.3f']; else ss=[ss sep '%g']; end
    doku=[doku dtab 'R'];
end
if isfield(N,'err')&&(length(N.err)==length(N.a))&&maxi,
    mess=[mess N.err(:)];
    ss=[ss sep '%g'];
    doku=[doku dtab 'err'];
end
if isfield(N,'konf')&&(length(N.konf)==length(N.a))&&maxi,
%     mess(1:length(N.k),end+1)=N.k(:);
    mess(1:length(N.konf),end+1)=N.konf(:);
    ss=[ss sep '%g'];
    doku=[doku dtab 'k'];%K/m
end
if isfield(N,'i')&&(length(N.i)==length(N.a))&&maxi,
    mess(1:length(N.i),end+1)=N.i(:);
    ss=[ss sep '%g'];%.4f
    doku=[doku dtab 'I'];%I/A
end
if isfield(N,'u')&&(length(N.u)==length(N.a))&&maxi,
    mess(1:length(N.u),end+1)=N.u(:);
    ss=[ss sep '%g'];%.4f
    doku=[doku dtab 'U'];%U/V
end
if isfield(N,'ip')&&(~isempty(find(N.ip)))&&(length(N.ip)==length(N.a))&&maxi,
    mess(1:length(N.ip),end+1)=N.ip(:);
    ss=[ss sep '%g'];%.3f
    doku=[doku dtab 'ip'];%IP/mrad
end
ss=[ss newline];
doku=[doku newline];
if maxi, fprintf(fid,doku); end
fprintf(fid,ss,mess');
if isfield(N,'topo')&&(length(N.topo)>0)&&maxi,
    fprintf(fid,['%d # Number of topographical points' newline],size(N.topo,1));
    fprintf(fid,['# x z positions for each topo point' newline]);
    fprintf(fid,['%g' sep '%g' newline],N.topo(:,1:2)');
end
fclose(fid);
