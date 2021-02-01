datfile='ro\ro-full.dat';
N=read3dfile(datfile);
[Model.x,Model.y,Model.z,Model.M,Model.Bg]=modelfromdata3d(N);
tic;S=calcsens3d(Model.x,Model.y,Model.z,N);toc
tic;Ss=spcalcsens3d(Model.x,Model.y,Model.z,N,1e-3);toc
tic;Ss1=spcalcsens3d1(Model.x,Model.y,Model.z,N,1e-3);toc
