function [R,Rrez]=collectrhoa(N,MEA,rhoBg)

% FDFWD3D - 3D DC Forward Calculation with finite differences
% Rhoa = collectrhoa(N,MEA,rhoBg)
% [Rhoa_abmn,Rhoa_mnab] = ...
%   N      - Structure of electrode numbers(a,b,m,n), 
%            k-factors(k) and measurements(r)
%            elec- Electrode Positions
%  MEA - multielectrode potential matrix
%        matrix of (size(N.elec,1))^2

if nargin<2, error('Two input arguments required!'); end
if nargin<3, rhoBg=0; end
% hacked!!!!!!!!!
if length(rhoBg)>1,
    for e=1:size(MEA,1),
        dx=N.elec(:,1)-N.elec(e,1);
        dy=N.elec(:,2)-N.elec(e,2);
        dxdy2=dx.^2+dy.^2;
        dz=N.elec(:,3)-N.elec(e,3);
        dz1=N.elec(:,3)+N.elec(e,3);
        MEA(:,e)=MEA(:,e)+(1./sqrt(dxdy2+dz.^2)+1./sqrt(dxdy2+dz1.^2))/4/pi/rhoBg(e);
%         ele(e,1)=1;%???
%         MEA(:,e)=MEA(:,e)+1./sqrt(sum(ele.^2,2))/2/pi*rhoBg(e);
    end
end

data=length(N.a);
R=zeros(size(N.a));Rrez=R;
for l = 1:data,
    R(l)=MEA(N.a(l),N.m(l));
    Rrez(l)=MEA(N.m(l),N.a(l));
    if N.n(l)>0, 
        R(l)=R(l)-MEA(N.a(l),N.n(l)); 
        Rrez(l)=Rrez(l)-MEA(N.n(l),N.a(l)); 
    end
    if N.b(l)>0,
        R(l)=R(l)-MEA(N.b(l),N.m(l));
        Rrez(l)=Rrez(l)-MEA(N.m(l),N.b(l));
        if N.n(l)>0, 
            R(l)=R(l)+MEA(N.b(l),N.n(l)); 
            Rrez(l)=Rrez(l)+MEA(N.n(l),N.b(l)); 
        end
    end
end
if isfield(N,'k'),
    R(:)=R(:).*N.k(:);
    Rrez(:)=Rrez(:).*N.k(:);
else
    display('warning! no geometric factors present in data file!');
end
if length(rhoBg)==1,
    R=R+rhoBg;
    Rrez=Rrez+rhoBg;
end
if nargout<2,
    rez=(R-Rrez)*2./(R+Rrez);
    messg(sprintf('Standard deviation of reciprocity %.2f%%',std(rez)*100));
    R=sqrt(abs(R.*Rrez));
end
