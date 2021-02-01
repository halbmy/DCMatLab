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
        ele=N.elec;
        ele(:,1)=ele(:,1)-N.elec(e,1);
        ele(:,2)=ele(:,2)-N.elec(e,2);
        ele(e,1)=1;
        MEA(:,e)=MEA(:,e)+1./sqrt(sum(ele.^2,2))/2/pi*rhoBg(e);
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
R(:)=R(:).*N.k(:);
Rrez(:)=Rrez(:).*N.k(:);
if length(rhoBg)==1,
    R=R+rhoBg;
    Rrez=Rrez+rhoBg;
end
if nargout<2,
    rez=(R-Rrez)*2./(R+Rrez);
    message(sprintf('Standard deviation of reciprocity %.2f%%',std(rez)*100));
    R=sqrt(abs(R.*Rrez));
end
