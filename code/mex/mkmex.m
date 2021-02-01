if ispc,
  display('Setup mex to Fortran compiler');
  mex -setup
end
outdir=['..' filesep '..' filesep 'dcapps' filesep 'dc2dinvres'];
mex('-O','-outdir',outdir,'sens2dxz.f');
outdir=['..' filesep '..' filesep 'dcapps' filesep 'dc3dinvres'];
mex('-O','-outdir',outdir,'sens3dfull.f');
mex('-O','-outdir',outdir,'sens3dplane.f');
if ispc,
  display('Setup mex to C compiler');
  mex -setup
end
mex('-outdir',outdir,'-inline','fd3dmea.c','ldl.c');
display('go to geo2dmx directory and run seperately and copy');
