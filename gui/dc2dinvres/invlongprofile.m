global N Mod INV MAL FIX FOR
xorg=Mod.x;FIX=[];
if min(diff(N.elec(:,1)))<0, error('sort electrodes before'); end
if ~isfield(INV,'rbnorm'), INV.rbnorm=1; end
if ~isfield(INV,'rbzfak'), INV.rbzfak=0.3; end
if ~isfield(INV,'lbound'), INV.lbound=0; end
if ~isfield(INV,'ubound'), INV.ubound=0; end

abmn=[N.a N.b N.m N.n];
ma=max(abmn,[],2);
mi=min(abmn,[],2);
dd=max(ma-mi)*2;
del=median(diff(N.elec(:,1)));

nel=size(N.elec,1);
rho=median(N.r);Mod.M(:)=rho;R=ones(size(N.r))*rho;
patch2dmodel(xorg,Mod.z,Mod.M,MAL,N);
l=0;mael=0;
while mael<nel,
    miel=fix(l*dd/2+1);
    mael=min(miel+dd-1,nel);if nel-mael<15, mael=nel; end
    NN.elec=N.elec(miel:mael,:);
    fi=find((mi>=miel)&(ma<=mael));
    NN.a=N.a(fi)-miel+1;NN.b=N.b(fi)-miel+1;NN.m=N.m(fi)-miel+1;NN.n=N.n(fi)-miel+1;
    NN.r=N.r(fi);NN.k=N.k(fi);NN.err=N.err(fi);
    ixa=max(max(find(xorg<NN.elec(1,1)))-2,1);
    ixe=min(min(find(xorg>NN.elec(end,1)))+2,length(xorg));
    fprintf('\nElectrodes %d-%d (%.1f..%.1fm) ',miel,mael,xorg(ixa),xorg(ixe));
    MM=Mod.M(ixa:ixe-1,:);Mod.x=xorg(ixa:ixe);nM=MM;RR=R(fi);
    if l>0,
        fixto=fix(size(MM,1)/4);fprintf('fixed to %.1fm ',Mod.x(fixto+1));
        FIX=zeros(size(MM));FIX(1:fixto,:)=-1;
    end
    fprintf('%.1f ',chi2(NN.r,RR,NN.err,INV.lolo));
    S=calcsens2d(Mod.x,Mod.z,NN);
    for i=1:3,
        dR=logtrans(NN.r,INV.lbound,INV.ubound)-logtrans(RR,INV.lbound,INV.ubound);
        dM0=logtrans(MM(:),INV.lbound,INV.ubound);
        dM=invers(S,dR,INV,NN.err,dM0);
        nM(:)=ilogtrans(logtrans(MM(:),INV.lbound,INV.ubound)+dM);
        R1=dcfwd2d(Mod.x,Mod.z,nM,Mod.Lay,NN,FOR);
        [fak,appR]=linesearch(NN,RR,R1,INV.lolo);
        if fak==0, fprintf('Line search failed (it. %d)!\n',i);break; end
        if fak<0.95,
            MM(:)=ilogtrans(logtrans(MM(:),INV.lbound,INV.ubound)+dM*fak);
            RR=dcfwd2d(Mod.x,Mod.z,MM,Mod.Lay,NN,FOR);
        else
            MM=nM;RR=R1;
        end
        chiq=chi2(NN.r,RR,NN.err,INV.lolo);
        fprintf('%.1f ',chiq);
        if chiq<1, break; end
    end
    Mod.M(ixa:ixe-1,:)=MM;
    patch2dmodel(xorg,Mod.z,Mod.M,MAL,N);pause(1.0);
    l=l+1;
end
Mod.x=xorg;

return
global N Mod.x Mod.z Mod.M INV MAL FIX Mod.Lay FOR
xorg=Mod.x;FIX=[];
if min(diff(N.elec(:,1)))<0, error('sort electrodes before'); end
if ~isfield(INV,'rbnorm'), INV.rbnorm=1; end
if ~isfield(INV,'rbzfak'), INV.rbzfak=0.3; end
if ~isfield(INV,'lbound'), INV.lbound=0; end
if ~isfield(INV,'ubound'), INV.ubound=0; end

abmn=[N.a N.b N.m N.n];
ma=max(abmn,[],2);
mi=min(abmn,[],2);
dd=max(ma-mi);
del=median(diff(N.elec(:,1)));

nel=size(N.elec,1);
rho=median(N.r);Mod.M(:)=rho;R=ones(size(N.r))*rho;
l=0;mael=0;dd=256;
while mael<nel,
    miel=l*dd/2+1;
    mael=min(miel+dd-1,nel);
    fprintf('\n%d %d ',miel,mael);
    NN.elec=N.elec(miel:mael,:);
    fi=find((mi>=miel)&(ma<=mael));
    NN.a=N.a(fi)-miel+1;NN.b=N.b(fi)-miel+1;NN.m=N.m(fi)-miel+1;NN.n=N.n(fi)-miel+1;
    NN.r=N.r(fi);NN.k=N.k(fi);NN.err=N.r(fi);
    ixa=max(max(find(xorg<NN.elec(1,1)))-2,1);
    ixe=min(min(find(xorg>NN.elec(end,1)))+2,length(xorg));
    MM=Mod.M(ixa:ixe-1,:);Mod.x=xorg(ixa:ixe);nM=MM;RR=R(fi);
    if l>0, FIX=zeros(size(MM));FIX(1:fix(size(MM,1)/4),:)=-1; end
    S=calcsens2d(Mod.x,Mod.z,NN);fprintf('%.1f ',chi2(NN.r,RR,NN.err,INV.lolo));
    for i=1:3,
        dR=logtrans(NN.r,INV.lbound,INV.ubound)-logtrans(RR,INV.lbound,INV.ubound);
        dM0=logtrans(MM(:),INV.lbound,INV.ubound);
        dM=invers(S,dR,INV,NN.err,dM0);
        nM(:)=ilogtrans(logtrans(MM(:),INV.lbound,INV.ubound)+dM);
        R1=dcfwd2d(Mod.x,Mod.z,nM,Mod.Lay,NN,FOR);
        [fak,appR]=linesearch(NN,RR,R1,INV.lolo);
        if fak==0, fprintf('Line search failed (it. %d)!\n',i);break; end
        if fak<0.95,
            MM(:)=ilogtrans(logtrans(MM(:),INV.lbound,INV.ubound)+dM*fak);
            RR=dcfwd2d(Mod.x,Mod.z,MM,Mod.Lay,NN,FOR);
        else
            MM=nM;RR=R1;
        end
        chiq=chi2(NN.r,RR,NN.err,INV.lolo);
        fprintf('%.1f ',chiq);
        if chiq<0.8, break; end
    end
    Mod.M(ixa:ixe-1,:)=MM;
    patch2dmodel(xorg,Mod.z,Mod.M,MAL,N);pause(1.0);
    l=l+1;
end
Mod.x=xorg;

return
A=zeros(size(N.elec,1));
for i=1:length(N.a),
    isb=(N.b(i)>0);isn=(N.n(i)>0);
    A(N.a(i),N.m(i))=1;
    if isb,
        A(N.b(i),N.m(i))=1;A(N.a(i),N.b(i))=1;
    end
    if isn,
        A(N.a(i),N.n(i))=1;
        if isb, A(N.b(i),N.n(i))=1; end
    end
    if isb&isn, A(N.b(i),N.n(i))=1; end
end
