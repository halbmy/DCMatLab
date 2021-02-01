function R = fefwdcall(N,resistivity,rhobg,dim)

% FEFWDCALL - call DCFEM (finite element forward routine)
% R = fefwdcall(N,resistivity[,rhobg,dim])
% (needs potentials in primaryPot/interpolated/potential & mesh/meshSec.bms)

if nargin<4, dim=2; end
if nargin<3, rhobg=median(resistivity); end
netz=['mesh' filesep 'meshSec.bms'];
potmat='primaryPot/interpolated/potential';
xtra='';
if dim==-2, 
    netz=['tmp' filesep 'meshSec.bms'];
    potmat=['tmp' filesep 'secPot.mat'];
    xtra='-B '; 
end % Circle geometry
rhomap=(1:length(resistivity)+1)';rhomap(1,2)=rhobg;
rhomap(2:end,2)=resistivity;save('rho.map','rhomap','-ascii')
fid=fopen('rho.map','w');fprintf(fid,'%d\t%e\n',rhomap');fclose(fid);
system(['dcfem -STHim -s2 -d' num2str(abs(dim)) ' -r' num2str(rhobg) ' -arho.map -x' potmat ' ' xtra netz]);
MEA=readcollect('pot.collect');R=collectrhoa(N,MEA);
