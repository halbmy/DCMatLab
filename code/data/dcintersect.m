function [Data,fi1,fi2,Data2a]=dcintersect(Data1,Data2,sortdata,delelecs)

% DCINTERSECT - Find intersecting data from 2 data sets
% [Data,index1,index2] = dcintersect(Data1,Data2)
% with Data.r = Data2.r/Data1.r OR
% [Data1a,Data2a] = dcintersect(Data1,Data2)
% [Data1a,index1,index2,Data2a] = dcintersect(Data1,Data2)

if nargin<4, delelecs=0; end
if nargin<3, sortdata=0; end
if sortdata,
    abmn1=[min(Data1.a(:),Data1.b(:)) max(Data1.a(:),Data1.b(:)) min(Data1.m(:),Data1.n(:)) min(Data1.m(:),Data1.n(:))];
    abmn2=[min(Data2.a(:),Data2.b(:)) max(Data2.a(:),Data2.b(:)) min(Data2.m(:),Data2.n(:)) min(Data2.m(:),Data2.n(:))];
else
    abmn1=[Data1.a(:) Data1.b(:) Data1.m(:) Data1.n(:)];
    abmn2=[Data2.a(:) Data2.b(:) Data2.m(:) Data2.n(:)];
end
[abmn,fi1,fi2]=intersect(abmn1,abmn2,'rows');
if delelecs,
    Data=deldeadelecs(extractmeasurement(Data1,fi1));
else
    Data=extractmeasurement(Data1,fi1);
end
if (nargout>3)|(nargout==2), % give back two data
    Data2a=deldeadelecs(extractmeasurement(Data2,fi2));
    if nargout==2, fi1=Data2a; end
else % give back the ratio of the 
    if isfield(Data1,'r')&&isfield(Data2,'r'),
        Data.r=Data2.r(fi2)./Data1.r(fi1);
    elseif isfield(Data1,'rho')&&isfield(Data2,'rho'),
        Data.r=Data2.rho(fi2)./Data1.rho(fi1);
    end
end