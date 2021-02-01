function NN=deldeadpos(N)

% DELDEADELECS - Delete dead electrodes
% Nnew = deldeadelecs(N)
% N, Nnew - data structure

map=zeros(size(N.pos,1),1);
used=unique([N.s;N.g]);
if used(1)==0, used(1)=[]; end
% notused=setxor(used,1:size(N.elec,1));
% N.elec(notused,:)=[];
NN=N;
NN.pos=N.pos(used,:);
if isfield(N,'x'), NN.x=N.x(used); end
if isfield(N,'y'), NN.y=N.y(used); end
if isfield(N,'z'), NN.z=N.z(used); end    
map(used)=1:length(used);
map=[0;map(:)];
NN.s=map(N.s+1);
NN.g=map(N.g+1);
if isfield(N,'t'), NN.t=N.t; end
