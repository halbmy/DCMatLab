function S = calcsens2d(x,z,N,INV)
% CALCSENS2D - Sensitivity calculation for 2D-cells
%              based on the fortran routine sens2dxz.dll
% S = calcsens2d(x,z,N)
% x,z .. Grid node vectors
% N .. struct electrode positions(elec),
%      numbers(a,b,m,n), rhoas(r) and k-factors(k)

t0=clock;
if (nargin>3)&&isfield(INV,'spsens')&&(INV.spsens>0),
    message('Sensitivity Integration (sparse) ...');    
    S=spcalcsens2dt(x,z,N,INV.spsens); 
    message(sprintf('ready(%ds) Elements: %s/%s(%.1f%%)',...
        round(etime(clock,t0)),int2hum(nnz(S)),int2hum(numel(S)),nnz(S)/numel(S)*100));
   return;
end
lx=length(x);
lz=length(z);

if (nargin>3)&&isfield(INV,'silent')&&(INV.silent>0), laut=0; else laut=1; end
%% Initialization
S=zeros(length(N.a),(lx-1)*(lz-1));
if laut, message('Sensitivity Integration...'); end
pause(0.1);
if laut, wb=waitbar(0,'Sensitivity integration'); end
data=length(N.r);
aller=fix(data/25);
mal=aller;
for l=1:data,
    sens=sens2dxz(x,z,N.elec(N.a(l),1:2),N.elec(N.m(l),1:2));
    if N.b(l)>0, sens=sens-sens2dxz(x,z,N.elec(N.m(l),1:2),N.elec(N.b(l),1:2)); end
    if N.n(l)>0, % A M N
        sens=sens-sens2dxz(x,z,N.elec(N.a(l),1:2),N.elec(N.n(l),1:2));
        if N.b(l)>0, sens=sens+sens2dxz(x,z,N.elec(N.n(l),1:2),N.elec(N.b(l),1:2)); end
    end
%     if N.n(l)>0, % A M N
%         sens=sens-sens2dxz(x,z,N.elec(N.a(l),1:2),N.elec(N.n(l),1:2));
%         if N.b(l)>0,
%             sens=sens-sens2dxz(x,z,N.elec(N.m(l),1:2),N.elec(N.b(l),1:2))+...
% 		 sens2dxz(x,z,N.elec(N.n(l),1:2),N.elec(N.b(l),1:2));
%         end
%     end
    S(l,:)=N.k(l)*sens(:)';
    mal=mal-1;
    if mal==0,
        if laut, waitbar(l/data); end
        mal=aller;
    end
end
if laut,
    close(wb);
    appendmessage(sprintf('ready(%.2f)',etime(clock,t0)));
end
