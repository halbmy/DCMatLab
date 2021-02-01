function N=readresecsfile(filename)

% READRESECSFILE - Read RESECS ascii export file (*.TXT)
% Data = readresecsfile(filename)
% Data..structure consisting of a,b,m,n,rhoa,U,I,...
% the field names will be interpreted from the header line

fid=fopen(filename,'r');
zeile=fgetl(fid);
if strcmp(lower(zeile(1:6)),'inject'),
    for i=1:4, zeile=fgetl(fid); end
end
i=0;irhoa=0;iip=0;ierr=0;itype=0;
ii=0;iu=0;ir=0;ik=0;
c1x=0;c1y=0;c1z=0;c2x=0;c2y=0;c2z=0;
p1x=0;p1y=0;p1z=0;p2x=0;p2y=0;p2z=0;
formstr='';%'%*s';
while ~isempty(zeile),
    i=i+1;
    [tok,zeile]=strtok(zeile);
    switch tok,
        case 'Type', formstr=[formstr '%*s'];i=i-1;
        case {'C1(x)','C1_X'}, c1x=i;
        case {'C1(y)','C1_Y'}, c1y=i;
        case {'C1(z)','C1_Z'}, c1z=i;
        case {'C2(x)','C2_X'}, c2x=i;
        case {'C2(y)','C2_Y'}, c2y=i;
        case {'C2(z)','C2_Z'}, c2z=i;
        case {'P1(x)','P1_X'}, p1x=i;
        case {'P1(y)','P1_Y'}, p1y=i;
        case {'P1(z)','P1_Z'}, p1z=i;
        case {'P2(x)','P2_X'}, p2x=i;
        case {'P2(y)','P2_Y'}, p2y=i;
        case {'P2(z)','P2_Z'}, p2z=i;
        case 'Rho', irho=i;
        case 'I', ii=i;
        case 'U', iu=i;
        case {'M','P'}, iip=i;
        case 'D', ierr=i;
    end
    if ~isempty(tok)
        if ~isequal(tok,'Type'), formstr=[formstr '%f']; end
    end
end
A=mytextscan(fid,formstr);
% zeile=fgetl(fid);
% while ~isempty(zeile),
%     while
%     zeile=fgetl(fid);
% end
fclose(fid);
% [c1x c1y c1z c2x c2y c2z p1x p1y p1z p2x p2y p2z]
if c1x==0, error('C1x missing'); end
AA=A{c1x};
if (c1y>0), AA=[AA A{c1y}]; else AA(:,2)=0; end
if c1z, AA(:,3)=A{c1z}; end
if c2x, 
    BB=A{c2x}; 
    if c2y, BB=[BB A{c2y}]; else BB(:,2)=0; end
else BB=[]; end
% if c2x*c2y>0, BB=[A{c2x} A{c2y}]; else BB=[]; end
if c2z, BB(:,3)=A{c1z}; end
if c1x==0, error('C1x missing'); end
MM=A{p1x};
if(p1y>0), MM=[MM A{p1y}]; else MM(:,2)=0; end
if p1z, MM(:,3)=A{c1z}; end
if p2x,
    NN=A{p2x};
    if p2y, NN=[NN A{p2y}]; else NN(:,2)=0; end
else NN=[]; end
% if p2x*p2y>0, NN=[A{p2x} A{p2y}]; else NN=[]; end
if p2z, NN(:,3)=A{c1z}; end
N.elec=unique(sortrows([AA;BB;MM;NN]),'rows');
data=length(A{1});
N.b=zeros(data,1);N.n=zeros(data,1);
[tf,N.a]=ismember(AA,N.elec,'rows');
if ~isempty(BB), [tf,N.b]=ismember(BB,N.elec,'rows'); end
[tf,N.m]=ismember(MM,N.elec,'rows');
if ~isempty(NN), [tf,N.n]=ismember(NN,N.elec,'rows'); end
if irho, N.r=A{irho}; end
if ii>0, N.i=A{ii}/1000; end
if iu>0, N.u=A{iu}/1000; end
if iip>0, N.ip=A{iip}; end
if ierr>0, N.err=A{ierr}*0.01; end
if ir>0, N.rho=A{ir}; else
    if ii*iu>0, 
        if any(N.i==0), 
            N=delmeasurement(N,N.i==0);
            display('deleted zero current data!');
        end        
        N.rho=N.u./N.i; 
    end
end
if (iu>0)&&any(N.u==0), 
    N=delmeasurement(N,N.u==0);
    display('deleted zero voltage data!');
end        
if iip&iu,
   fi=find(N.u<0);
   N.ip(fi)=pi*1000-N.ip(fi);
end
if ik>0, N.k=A{ik}; else 
    if size(N.elec,2)>2, N.k=getkonf3d(N); else N.k=getkonf2d(N); end
end
if irhoa>0, N.r=A{irhoa}; else
    if isfield(N,'rho')&isfield(N,'k'), N.r=N.rho.*N.k; end
end
% if any(N.r<=0), 
%     N=delmeasurement(N,N.r<=0);
%     display('deleted negative/zero apparent resistivities!');
% end        