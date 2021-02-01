function [X,Y,Z,M,Bg]=modelfromdata3d(N,dx,dy)
% MODELFROMDATA - Create 2d-model from data
% Model = modelfromdata(N[,dx[,dy]])
% [x,y,z,M,Bg] = modelfromdata(N[,dx[,dy]])
% N..Structure of electrode numbers(a,b,m,n), 
%    k-factors(k) and measurements(r)
%    elec..Electrode Positions
% dx/dy..Grid Spacings
% x,y,z..Grid nodes
% M..Model resistivities
% Bg..Background resistivities

message(sprintf('Creating grid model from Data'));
if (nargin<3)||(dy==0),
    aa=diff(unique(sort(N.elec(:,2))));
    dy=min(aa(find(aa)));dy=round(dy*1000)/1000;
    if dy<0.5, 
        if isfield(N,'zweid'),
            dy=min(diff(unique(N.zweid{1}.elec(:,1))));    
        else
            dy=1; 
        end
    end
    if dy>0.3, dy=round(dy*20)/20; else dy=rndig(dy,2); end
end
if (nargin<2)||(dx==0),
    aa=diff(unique(sort(N.elec(:,1))));
    dx=min(aa(find(aa)));
    dx=round(dx*1000)/1000;
    if dx<0.5, 
        if isfield(N,'zweid'),
            dx=min(diff(unique(N.zweid{1}.elec(:,1))));    
        else
            dx=1; 
        end
    end
    if dx>0.3, dx=round(dx*20)/20; else dx=rndig(dx,2); end
else
    dy=dx;
end
if isempty(dx), dx=dy; end
if isempty(dy), dy=dx; end

minx=floor(min(N.elec(:,1))/dx)*dx;
maxx=floor(max(N.elec(:,1))/dx+1)*dx;
miny=floor(min(N.elec(:,2))/dy)*dy;
maxy=floor(max(N.elec(:,2))/dy+1)*dx;

%message(sprintf('Min(x)=%g Max(x)=%g Min(y)=%g Max(y)=%g D=%g',minx,maxx,miny,maxy,D)); 
zusatz=1;
X=minx-zusatz*dx:dx:maxx+zusatz*dx;
Y=miny-zusatz*dy:dy:maxy+zusatz*dy;
if isfield(N,'d'),
    dz=sqrt(dx*dy);
    Z=0:dz:max(abs(N.d))+dz;
else    
    Z=zparam(N);
end
X=round(X*1000)/1000;Y=round(Y*1000)/1000;
I=length(X);J=length(Y);K=length(Z);
if I>101, X=linspace(min(X),max(X),51);I=length(X); end
if J>101, Y=linspace(min(Y),max(Y),51);J=length(Y); end
message(sprintf('I=%d J=%d K=%d => %d Cells',I-1,J-1,K-1,...
    (I-1)*(J-1)*(K-1)));
message(sprintf('X=%.1f-%.1f(%.1f) Y=%.1f-%.1f(%.1f) Z=(%.1f-%.1f)',X(1),X(end),dx,Y(1),Y(end),dy,Z(1),Z(end)));

% Sigma_q=Mittelwert der geringsten Eindringtiefe o. Mittel
% minKonf=min(abs(N.k)); 
% rq=median(N.r(find(abs(N.k)==minKonf)));
rq=median(N.r); %!!!
if rq>20, rq=round(rq/10)*10; end
message(sprintf('Rho_1 set to %.1f',rq));
M=ones(I-1,J-1,K-1)*rq;
Bg=ones(K,1)*rq;
R=ones(size(N.a))*rq;
if nargout==1,
    X=struct('x',X,'y',Y,'z',Z,'M',M,'Bg',Bg,'R',R); 
end