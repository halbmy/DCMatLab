function NN=deldeadelecs(N)

% DELDEADELECS - Delete dead electrodes
% Nnew = deldeadelecs(N)
% N, Nnew - data structure

map=zeros(size(N.elec,1),1);
used=unique([N.a;N.b;N.m;N.n]);
if used(1)==0, used(1)=[]; end
% notused=setxor(used,1:size(N.elec,1));
% N.elec(notused,:)=[];
NN=N;
NN.elec=N.elec(used,:);
if isfield(N,'x'), NN.x=N.x(used); end
if isfield(N,'y'), NN.y=N.y(used); end
if isfield(N,'z'), NN.z=N.z(used); end    
map(used)=1:length(used);
map=[0;map(:)];
NN.a=map(N.a+1);
NN.b=map(N.b+1);
NN.m=map(N.m+1);
NN.n=map(N.n+1);
if isfield(N,'r'), NN.r=N.r; end
if isfield(N,'k'), NN.k=N.k; end
if isfield(N,'err'), NN.err=N.err; end
if isfield(N,'ip'), NN.ip=N.ip; end
if isfield(N,'i'), NN.i=N.i; end
if isfield(N,'u'), NN.u=N.u; end
if isfield(N,'rez'), NN.rez=N.rez; end
if isfield(N,'rho'), NN.rho=N.rho; end
if isfield(N,'zweid'), NN.zweid=N.zweid; end
if isfield(N,'nr'), NN.nr=N.nr; end
if isfield(N,'names'), NN.names=N.names; end