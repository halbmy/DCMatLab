function xz=mbm2xz(xmb,topo,topoisxz,rightpoint)

% MBM2XZ - Converts tape measure coordinates in topographical xz
% xz = mbm2xz(xtape,topography[,topoisxz,rightpoint])
% xz         .. positions in real coordinates
% xmb        .. tape measure x coordinate
% topo       .. topography (paired values)
% topoisxz   .. topo is x z (otherwise tape-x z) [default=0]
% rightpoint .. number of point whose xmb equals topox [default=1]
% Thomas Günther (thomas@resistivity.net)

if nargin<4, rightpoint=1; end
if nargin<3, topoisxz=0; end

if topoisxz, % convert first topography_x to mbm
    topox=topo(1,1)+[0;cumsum(sqrt(sum(diff(topo).^2,2)))];
    shift=topo(rightpoint,1)-topox(rightpoint);
    topo(:,1)=topox+shift;
else
    shift=0;
end
z=interp1(topo(:,1),topo(:,2),xmb,'linear','extrap');
x=xmb(1)+[0;cumsum(sqrt(diff(xmb).^2-diff(z).^2))];
xz=[x(:)-shift z(:)];
di=sqrt(sum(diff(xz).^2,2))';
fprintf('min/max electrode distance are %.2f/%.2fm\n',min(di),max(di));