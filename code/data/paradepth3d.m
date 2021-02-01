function zpar=zparam(N,frac)

% PARDEPTH3D - Computes penetration depth by 1d sensitivity
% z_par = paradepth3d(N);
% N=struct of arrays a,b,m,n(electrode numbers), r(rhoa), k(konf)
%             elec..Coordinates of Electrodes

if nargin<2, frac=0.95; end
amq=[];bmq=amq;anq=amq;bnq=amq;
data=length(N.a);
fb=find(N.b);fn=find(N.n);fbn=find(N.b.*N.n);
amq=sum((N.elec(N.a,:)-N.elec(N.m,:)).^2,2);
anq=amq;bmq=amq;bnq=amq;
bmq(fb)=sum((N.elec(N.b(fb),:)-N.elec(N.m(fb),:)).^2,2);
anq(fn)=sum((N.elec(N.a(fn),:)-N.elec(N.n(fn),:)).^2,2);
bnq(fbn)=sum((N.elec(N.b(fbn),:)-N.elec(N.n(fbn),:)).^2,2);
% wieso ist das auskommentiert?
% amq=0.25*amq;bmq=0.25*bmq;anq=0.25*anq;bnq=0.25*bnq;
bnq(find(isnan(bnq)))=Inf;
am=sqrt(amq);an=sqrt(anq);bm=sqrt(bmq);bn=sqrt(bnq);
depth=0.1;
ss=0;
zi=[];
si=[];
sens=0;
scal=1./am;scal(fn)=scal(fn)-1./an(fn);
scal(fb)=scal(fb)-1./bm(fb);scal(fbn)=scal(fbn)+1./bn(fbn);
while sens<0.99,
    %zzq=repmat(4*depth.^2,data,1);
    zq=4*depth^2;
    ss=1./sqrt(amq+zq);
    ss(fn)=ss(fn)-1./sqrt(anq(fn)+zq);
    ss(fb)=ss(fb)-1./sqrt(bmq(fb)+zq);
    ss(fbn)=ss(fbn)+1./sqrt(bnq(fbn)+zq);
    %sens=sqrt(sum(((1-ss./scal).^2)/data));
    sens=sum(1-ss./scal)/data;
    %ss=am./sqrt(amq+zzq)-an./sqrt(anq+zzq)-bm./sqrt(bmq+zzq)+bn./sqrt(bnq+zzq);
    %sens=sqrt(sum((ss).^2)/data);
    zi=[zi depth];
    si=[si sens];
    depth=depth*1.1;
end
[mi,nb]=min(abs(si-frac));
zpar=zi(nb);
if nargout==0, 
    plot(zi,si,'bx-'); 
    line(xlim,frac*[1 1],'Color','red');
    line(zpar(end)*[1 1],ylim,'Color','red');
    title(sprintf('%.1f%% sensitivity at %.2fm',frac*100,zpar(end)));
end