function [x,z,M]=modelfromdata2d(N,zusatz)

% MODELFROMDATA - Create 2d-model from data
% [x,z,M] = modelfromdata(N[,OPT])
% N..Structure of Electrode numbers(a,b,m,n), 
% k-factors(k), measurements(r) and
% N.elec..Electrode Positions
% x/z vectors of grid nodes
% M - Model matrix ( length(x)-1 * length(z) )

if nargin<2, zusatz=1; end
del=diff(sort(N.elec(:,1)));
% dx=min(del(find(del>1e-3)));
dx=median(del(del>1e-3));
dx=round(dx*1000)/1000; %rounding on mm
% dx=1;zusatz=1;
%dx=dx/2; % only a test!!!
minx=min(N.elec(:,1));
maxx=max(N.elec(:,1));
if isempty(dx), 
    dx=(maxx-minx)/25; 
    l=1;while dx<1, dx=dx*10;l=l*10; end
    dx=round(dx)/l;
end
if isnan(dx), dx=1; end
x=minx-dx*zusatz:dx:maxx+dx*zusatz;
zel=N.elec(:,2);
mzel=max(zel);
if mzel>0,
    dzel=min(diff(unique(zel)));
    if (isempty(dzel))||(length(unique(N.elec(:,1)))>size(N.elec,1)/2), % underwater measurement
        z=zparam(N); % from 1D-Sens.
        zm=mean(zel);
        zunten=zparam(N,length(z),max(z));
        if isempty(dzel), % all constant
            zoben=linspace(0,zm,round(zm/dx)+1);
        else % underwater topography
            zoben=linspace(0,zm,round(zm/zunten(2))+1);
        end
        z=[zoben(1:end-1) zunten+zm];
    else % borehole2surface
        z=0:dzel:mzel+dzel*zusatz;
    end
else
    z=zparam(N); % from 1D-Sens.
    z=zparam(N,length(z),max(z));
end
%z=2*z; % only a test
%zneu=round(10*(z(end)*z(end)/z(end-1)))/10;
%z=[z zneu];
% Prolongation
x=round(x(:)*1000)/1000;
z=round(z(:)*1000)/1000;
%minkonf=min(abs(N.k(:)));
%rq=median(N.r(abs(N.k)==minkonf));
rq=median(N.r);
if rq>20, rq=round(rq/10)*10; end
% rq=100;
%M=ones(length(z)-1,length(x)-1)*rq;
M=ones(length(x)-1,length(z)-1)*rq;
messg(sprintf('Creating new model (%d cells) Rho_0 = %.1f',numel(M),rq));
messg(strcat(sprintf('x=%.1f..(%.1f)..%.1f',min(x),dx,max(x)),...
    '  z= ',sprintf('%.1f ',z)));
if nargout==1,
    xx=x;
    x=[];
    x.x=xx;
    x.z=z;
    x.M=M;
    x.Lay=rq;
    x.R=ones(size(N.r))*rq;
    x.Cov=[];
end