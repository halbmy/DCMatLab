Model.x=0:3;
Model.y=-1:1;
Model.z=0:4;
Model.M=ones(length(Model.x)-1,length(Model.y)-1,length(Model.z)-1);
Data=[];Data.elec=[1 0 1;1.5 0 2;2 0 1;2 0 2];
Data.a=1;Data.b=2;Data.m=3;Data.n=4;
Data.k=getkonf3d(Data);
FOR=struct('method',0,'direct',1);
rhoa=mfdfwd3d(Model,Data,FOR)
return
SIGMA=ones(length(Model.x)+1,length(Model.y)+1,length(Model.z)+1);
SIGMA(2:end-1,2:end-1,2:end-1)=1./Model.M;SIGMA(:,:,end)=SIGMA(:,:,end-1);
SIGMA(1,:,:)=SIGMA(2,:,:);SIGMA(end,:,:)=SIGMA(end-1,:,:);
SIGMA(:,1,:)=SIGMA(:,2,:);SIGMA(:,end,:)=SIGMA(:,end-1,:);
C=diskr(Model.x,Model.y,Model.z,SIGMA);C1=diskr(Model.x,Model.y,Model.z,(SIGMA>0));
PM=potmap(Data.elec,Model.x,Model.y,Model.z);
Rbg=ones(size(Data.elec,1),1);p=symamd(C);
MEA=fd3dmea(Model.x,Model.y,Model.z,Data.elec,C,C1,PM,Rbg,p)

