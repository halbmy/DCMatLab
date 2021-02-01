if ispc,
  display('Setup mex to Fortran compiler');
  mex -setup
end
<<<<<<< mkmex.m
outdir=fullfile('..','..','dcapps','dc2dinvres');
=======
outdir=['..' filesep '..' filesep 'dcapps' filesep 'dc2dinvres'];
>>>>>>> 1.3
mex('-O','-outdir',outdir,'sens2dxz.f');
<<<<<<< mkmex.m
outdir=fullfile('..','..','dcapps','dc3dinvres');
=======
outdir=['..' filesep '..' filesep 'dcapps' filesep 'dc3dinvres'];
>>>>>>> 1.3
mex('-O','-outdir',outdir,'sens3dfull.f');
mex('-O','-outdir',outdir,'sens3dplane.f');
if ispc,
  display('Setup mex to C compiler');
  mex -setup
end
mex('-outdir',outdir,'-inline','fd3dmea.c','ldl.c');
display('go to geo2dmx directory and run seperately and copy');
cd('geo2d');
mkgeo2dmx;
cd('..');
display('ready with mexing files');
