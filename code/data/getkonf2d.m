function konf=getkonf2d(N)

% GETKONF - Get Configuration factor from electrode positions
% N.k = getkonf(N)
% N.....structure of * arrays a,b,m,n = electrode numbers(elec)
%                    * elec = Electrode positions ( x,z )
% needed to set element k = konfiguration factor for N

ael=ones(length(N.a),size(N.elec,2))*Inf;
mel=-ones(length(N.a),size(N.elec,2))*Inf;
fi=find(N.a);ael(fi,:)=N.elec(N.a(fi),:);
fi=find(N.m);mel(fi,:)=N.elec(N.m(fi),:);
ast=ael;ast(:,2)=-ast(:,2);
% nn=find(N.a);
% if ~isempty(nn),
    konf=1./sqrt(sum((ael-mel).^2,2))+...
        1./sqrt(sum((ast-mel).^2,2));
% end
nn=find(N.n~=0);
if ~isempty(nn),
    nel=N.elec(N.n(nn),:);
    konf(nn)=konf(nn)-1./sqrt(sum((ael(nn,:)-nel).^2,2))-...
        1./sqrt(sum((ast(nn,:)-nel).^2,2));
end
nn=find(N.b~=0);
if ~isempty(nn),
    bel=N.elec(N.b(nn),:);
    bst=bel;bst(:,2)=-bst(:,2);
    konf(nn)=konf(nn)-1./sqrt(sum((bel-mel(nn,:)).^2,2))-...
    1./sqrt(sum((bst-mel(nn,:)).^2,2));
end
nn=find(N.b.*N.n~=0);
if ~isempty(nn),
    bel=N.elec(N.b(nn),:);
    bst=bel;bst(:,2)=-bst(:,2);
    nel=N.elec(N.n(nn),:);
    konf(nn)=konf(nn)+1./sqrt(sum((bel-nel).^2,2))+...
        1./sqrt(sum((bst-nel).^2,2));
end
ep=1e-12;
konf=round(4*pi./konf/ep)*ep;
