function N=delmeasurement(N,fi)

% DELMEASUREMENT - Delete measurement from data struct
% N = delmeasurement(N,fi)

if nargin<2, return; end
if min(fi)<=0, fi=find(fi); end
fie=fieldnames(N);
if isfield(N,'a'), ndata=length(N.a); 
elseif isfield(N,'s'), ndata=length(N.s); 
else fprintf('Cannot detect valid data');return; end
for i=1:length(fie),
   ff=getfield(N,fie{i}); 
   if (min(size(ff))==1)&&(length(ff)==ndata),
      ff(fi)=[];
      N=setfield(N,fie{i},ff);
   end
end