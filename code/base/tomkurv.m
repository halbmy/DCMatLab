function ind=tomkurv(rho,eta,ak)

% TOMKURV - L-curve criterion by maximum curvature
% index = tomkurv(rho,eta,lambdas)

if nargin<3, load rhoeta; end
%cc=kruemm(rho,eta,sqrt(ak));
cc=-curvature(rho,eta,0);
fak=max(eta)/max(abs(cc));
dc=diff(cc);
pr=-dc(1:end-1).*dc(2:end);
%ind=find(pr==max(pr))+1;
ind=find(pr>0)+1;
if isempty(ind),
    ind=1;
else
    [pp,j]=max(cc(ind));
    j=1; % insurance
    ind=ind(j);
end
if nargout==0,
    clf
    plot(rho,eta,'bx');
    hold on
    plot(rho,cc*fak,'go');
    plot(rho(ind),eta(ind),'r+');
    plot(rho(ind),cc(ind)*fak,'r+');
    hold off
end