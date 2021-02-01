function S=sens1d(Z,N)

% SENS1D - Sensitivity of 1D parameterization
% S = sens1d(z,N);

data=length(N.r);
amq=zeros(data,1);bmq=amq;anq=amq;bnq=amq;
for l = 1:data,
    aaa=N.elec(N.a(l),:);
    mmm=N.elec(N.m(l),:);
    if N.b(l)>0, bbb=N.elec(N.b(l),:); else bbb=[Inf Inf 0]; end
    if N.n(l)>0, nnn=N.elec(N.n(l),:); else nnn=[Inf Inf 0]; end
    amq(l)=sum((aaa-mmm).^2);    
    anq(l)=sum((aaa-nnn).^2);    
    bmq(l)=sum((bbb-mmm).^2);    
    bnq(l)=sum((bbb-nnn).^2);    
end
%amq=0.25*amq;bmq=0.25*bmq;anq=0.25*anq;bnq=0.25*bnq;
bnq(find(isnan(bnq)))=Inf;

scal=1./sqrt(amq)-1./sqrt(anq)-1./sqrt(bmq)+1./sqrt(bnq);
S=zeros(length(N.a),length(Z));
ss=scal;
for kk = 2:length(Z),
  zzq=ones(data,1)*(4*Z(kk)^2);
  oss=ss;
  ss=1./sqrt(amq+zzq)-1./sqrt(anq+zzq)-1./sqrt(bmq+zzq)+1./sqrt(bnq+zzq);
  sens=(oss-ss)./scal; % obere Grenze-untere Grenze
  S(:,kk-1)=sens;
end
S(:,end)=ss./scal;
return
%if INV.lolo==2, soll=soll*Bg(kk+1)./N.r';  %% Korrektur für d ln(rho)
%end
% if INV.lolo==2
%   suso=sum(Soll,2);
%   for kk = 1:size(Soll,2),
%      Soll(:,kk)=Soll(:,kk)./suso(:);
%   end
% end