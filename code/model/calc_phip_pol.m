function phip=calc_phip(X,Y,Z,sq,RR,RA)

IA=1-0.5*(RA(3)>0); % Source on/below surface

I=length(X);J=length(Y);K=length(Z);
IJK=I*J*K;

if exist('RR')~=1,
    RR=[reshape(repmat(X,1,J*K),1,IJK);reshape(repmat(repmat(Y,1,K),I,1),1,IJK);reshape(repmat(Z,I*J,1),1,IJK)];
end
Ra=sqrt(sum((RR-repmat(reshape(RA,3,1),1,IJK)).^2));
fi=find(Ra==0);
Ra(fi)=sqrt(min(diff(X))*min(diff(Y)))/2;
phip=IA./Ra/(2*pi*sq);
% phip=zeros(size(Ra));
% fa=find(Ra);
% phip(fa)=IA./Ra(fa)/(2*pi*sq);

% When electrode on node!
werte=[];
for ll = find(abs(phip)==Inf),
    [i,j,k]=antiindex(ll,I,J);
    px=phip(ll+1);xp=phip(ll-1);
    py=phip(ll+I);yp=phip(ll-I);
    plus=0;minus=0;lll=1;
    if isfinite(px*xp),
        plus=plus+1/(X(i+1)-X(i))+1/(X(i)-X(i-1));
        minus=minus+px+xp;
        lll=lll*2;
    end
    if isfinite(py*yp),
        plus=plus+1/(Y(j+1)-Y(j))+1/(Y(j)-Y(j-1));
        minus=minus+py+yp;
        lll=lll*2;
    end
    werte=[werte (sign(phip(ll))*plus-minus)/lll];
end
werte=werte*2/(pi*sq);
lll=0;
for ll = find(abs(phip)==Inf),
    lll=lll+1;
    phip(ll)=werte(lll);
end

phip=reshape(phip,I,J,K);
