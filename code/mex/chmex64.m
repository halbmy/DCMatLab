mex -o sens2dxz sens2dxz64.f
a=sens2dxz(0:1,0:1,[0 0],[1 0])
if a+0.014<1e-4, display('sens2dxz passed.'); end
mex sens3dfull.f
a=sens3dfull(0:1,0:1,0:1,[0 0 0],[1 0 0])
if a+0.011<1e-4, display('sens3dfull passed.'); end
%%
mex -largeArrayDims -inline -output fd3dmea fd3dmea64.c ldl.c
Data=[];Data.elec=[0 0 0;1 0 0];
Data.a=1;Data.m=2;Data.b=0;Data.n=0;
Data.k=getkonf3d(Data);
Model=[];Model.x=-3:4;Model.y=Model.x;Model.z=0:4;
Model.M=ones(length(Model.x)-1,length(Model.y)-1,length(Model.z)-1)*200;
Model.M(:,:,1)=100;Model.M(1,1,2)=201;
Model.Bg=ones(size(Model.z))*200;Model.Bg(1)=100;
R=mfdfwd3d(Model,Data)
if (R>100)&(R<200), display('fd3dmea passed.'); end
