function [N,index]=combdata2d(N,N1)

% COMBDATA2D - Combine data sets
% N = combdata2d(N1,N2)
% N,N1,N2..structures of * arrays a,b,m,n = electrode numbers(elec)
%                                 k = konfiguration factor
%                        * elec..Electrode positions ( x,z )


data=length(N.a);
ne=size(N.elec,1);
% ne1=size(N.elec,1);
data1=length(N1.a);

index=(1:size(N1.elec,1))'+ne;
%Elektroden anhängen
[aa,bb]=meshgrid(N1.elec(:,1)+N1.elec(:,2)*12.34,N.elec(:,1)+N.elec(:,2)*12.34);
[ii,jj]=find((aa-bb)==0);
index(jj)=ii;
ind=find(index>ne);
N.elec=[N.elec;N1.elec(ind,:)];
cu=cumsum(index>ne);
index(ind)=ne+cu(ind);
N.a(data+1:data+data1)=index(N1.a);
N.b(data+data1)=0; % verlängern
fb=find(N1.b>0);
N.b(fb+data)=index(N1.b(fb));
N.m(data+1:data+data1)=index(N1.m);
N.n(data+data1)=0; % verlängern
fn=find(N1.n>0);
N.n(fn+data)=index(N1.n(fn));
if isfield(N,'r')&isfield(N1,'r'),
    N.r=[N.r(:);N1.r(:)];
else N.r=[]; end
if isfield(N,'k')&isfield(N1,'k'),
    N.k=[N.k(:);N1.k(:)];
else N.k=[]; end
if isfield(N,'err')&isfield(N1,'err'),
    N.err=[N.err(:);N1.err(:)];
else N.err=[]; end
if isfield(N,'rho')&isfield(N1,'rho'),
    N.rho=[N.rho(:);N1.rho(:)];
else N.rho=[]; end
if isfield(N,'ip')&isfield(N1,'ip'),
    N.ip=[N.ip(:);N1.ip(:)];
else N.ip=[]; end
if isfield(N,'i')&isfield(N1,'i'),
    N.i=[N.i(:);N1.i(:)];
else N.i=[]; end
if isfield(N,'u')&isfield(N1,'u'),
    N.u=[N.u(:);N1.u(:)];
else N.u=[]; end
N=sort2delecs(N,1);
% index
% N
% size(N.elec)
% allei=1:length(N.a);
% plot(N.elec(N.a,1),allei,'r.');
% hold on
% plot(N.elec(N.m,1),allei,'b.');
% fb=find(N.b>0);
% plot(N.elec(N.b(fb),1),allei(fb),'r.');
% fn=find(N.n>0);
% plot(N.elec(N.n(fn),1),allei(fn),'b.');
% hold off