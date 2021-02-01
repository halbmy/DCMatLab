global Pro progdir
if isfield(Pro,'dirname'), cd(Pro.dirname); end
S=loadsens;
mm=load('smoothness.matrix');
C=sparse(mm(:,1)+1,mm(:,2)+1,mm(:,3));
D=spdiags(1./log(N.err+1),0,length(N.err),length(N.err));
NN=readcirc2dfile('bam.dat');
ch=cglscdp(S,NN.ip,Pro.lambda,C,D);