function N=extractmeasurement(N,fi)

% EXTRACTMEASUREMENT - Extract measurement from data struct
% NewData = extractmeasurement(Data,fi)
% Data .. data struct with posisions and fields
% fi .. vector of indices or logical values
% example: Data1 = extractmeasurement(Data,Data.u>0)
% See also: delmeasurement

if nargin<2, return; end
if min(fi)<=0, fi=find(fi); end
fie=fieldnames(N);
if isfield(N,'a'),
    ndata=length(N.a);
elseif isfield(N,'s'),
    ndata=length(N.s); 
else
    error('did not found a or s field!');
end
for i=1:length(fie),
   ff=getfield(N,fie{i}); 
   if (min(size(ff))==1)&&(length(ff)==ndata),
       N=setfield(N,fie{i},ff(fi));
   end
end
