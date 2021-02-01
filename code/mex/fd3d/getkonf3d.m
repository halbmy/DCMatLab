function konf=getkonf3d(N,negz)

% GETKONF - Get Configuration factor from electrode positions
% N.k = getkonf(N)
% N.....structure of * arrays a,b,m,n = electrode numbers(elec)
%                                  k = konfiguration factor
%                    * elec = Electrode positions ( x,z )

if nargin<2, negz=0; end
if negz, N.elec(:,3)=0; end
ael=N.elec(N.a,:);mel=N.elec(N.m,:);
ast=ael;ast(:,3)=-ast(:,3);
konf=1./sqrt(sum((ael-mel).^2,2))+...
    1./sqrt(sum((ast-mel).^2,2));
nn=find(N.n~=0);
if ~isempty(nn),
    nel=N.elec(N.n(nn),:);
    konf(nn)=konf(nn)-1./sqrt(sum((ael(nn,:)-nel).^2,2))-...
        1./sqrt(sum((ast(nn,:)-nel).^2,2));
end
nn=find(N.b~=0);
if ~isempty(nn),
    bel=N.elec(N.b(nn),:);
    bst=bel;bst(:,3)=-bst(:,3);
    konf(nn)=konf(nn)-1./sqrt(sum((bel-mel(nn,:)).^2,2))-...
    1./sqrt(sum((bst-mel(nn,:)).^2,2));
end
nn=find(N.b.*N.n~=0);
if ~isempty(nn),
    bel=N.elec(N.b(nn),:);
    bst=bel;bst(:,3)=-bst(:,3);
    nel=N.elec(N.n(nn),:);
    konf(nn)=konf(nn)+1./sqrt(sum((bel-nel).^2,2))+...
        1./sqrt(sum((bst-nel).^2,2));
end
ep=1e-12;
konf=round(4*pi./konf/ep)*ep;
