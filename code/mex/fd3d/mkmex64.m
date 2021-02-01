%if ispc,
  %display('Setup mex to Fortran compiler');
  %mex -setup
%end
%outdir=fullfile('..','..','dcapps','dc2dinvres');
%mex('-O','-outdir',outdir,'sens2dxz.f');
outdir=fullfile('..','..','dcapps','dc3dinvres');
%mex('-O','-outdir',outdir,'sens3dfull.f');
%mex('-O','-outdir',outdir,'sens3dplane.f');
%if ispc,
  %display('Setup mex to C compiler');
  %mex -setup
%end
mex('-largeArrayDims','-DLDL_LONG','-outdir',outdir,'-output','fd3dmea','-inline','fd3dmea64.c','ldl64.c');
%display('go to geo2dmx directory and run seperately and copy');
%cd('geo2d');
%mkgeo2dmx;
%cd('..');
%display('ready with mexing files');
