outdir='.';
mex('-largeArrayDims','-DLDL_LONG','-outdir',outdir,'-output','fd3dmea','-inline','fd3dmea64.c','ldl64.c');
addpath('/home/guenther/src/dcmatlab/base');
addpath('/home/guenther/src/dcmatlab/model');
addpath('/home/guenther/src/dcmatlab/data');
N.elec=[0 0;1 0;0 1;1 1];N.elec(:,3)=0;
N.a=[1;1];N.b=[2;3];N.m=[3;2];N.n=[4;4];
N.k=getkonf3d(N);
Mod.x=-2:3;
Mod.y=-2:3;
Mod.z=0:3;
Mod.M=ones(length(Mod.x)-1,length(Mod.y)-1,length(Mod.z)-1)*100;
Mod.M(:,:,3:end)=200;
Mod.Bg=0;
R=mfdfwd3d(Mod,N)
