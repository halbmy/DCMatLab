function pot=potxy(x,y,sq,Phis,X,Y,RA,RB)
%% POTXY - Compute total potential from normal and anomal parts
%% pot = potxy(x,y,Phis,X,Y,RA,RB) where x and y are vectors of coordinates
%% needs RA RB X Y Phis sq as gloabl variables

I=length(X);
J=length(Y);
IA=1-0.5*(RA(3)>0); % Source on/below surface
IB=0;
if nargin>7, IB=-1+0.5*(RB(3)>0);  end

RR=[x(:) y(:)]';
RR(3,1)=0;
lx=length(x);
ra=sqrt(sum((RR-repmat(RA(:),1,lx)).^2));
if IB~=0, rb=sqrt(sum((RR-repmat(RB(:),1,lx)).^2)); end
fa=find(ra);
phip=zeros(size(ra));
phip(fa)=IA./ra(fa);
if IB~=0,
  fb=find(rb);
  phip(fb)=phip(fb)+IB./rb(fb);
end
while length(phip)<lx, phip=[phip 0]; end
pot=phip/(2*pi*sq);
k=1;
for nn=1:length(x),
  %i=1;j=1;k=1;
  %while (X(i)<x(nn))&(i<I), i=i+1; end
  %while (Y(j)<y(nn))&(j<J), j=j+1; end
  i=max(find(X<x(nn)));
  j=max(find(Y<y(nn)));
%  if (i<I)&(j<J),
      phis=Phis(i+1,j+1,k)*(x(nn)-X(i))*(y(nn)-Y(j))+Phis(i,j+1,k)*(X(i+1)-x(nn))*(y(nn)-Y(j))+Phis(i+1,j,k)*(x(nn)-X(i))*(Y(j+1)-y(nn))+Phis(i,j,k)*(X(i+1)-x(nn))*(Y(j+1)-y(nn));
      phis=phis/(X(i+1)-X(i))/(Y(j+1)-Y(j));
      pot(nn)=pot(nn)+phis;
      %  else
%      fprintf('Index too big! i=%d j=%d\n',i,j);
%  end
end
%return;
%fprintf('(%g;%g) found in (%g..%g);(%g..%g)\n',x,y,X(i),X(i+1),Y(j),Y(j+1));
