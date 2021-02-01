function Nnew=getmeasurement(N,fi)

% GETMEASUREMENT - Delete measurement from data struct
% N = getmeasurement(N,fi)

if nargin<2, return; end
if min(fi)<=0, fi=find(fi); end
fie=fieldnames(N);
ndata=length(N.a);
Nnew=[];
for i=1:length(fie),
   ff=getfield(N,fie{i}); 
   if (min(size(ff))==1)&&(length(ff)==ndata),
      Nnew=setfield(Nnew,fie{i},ff(fi)');
   end
end