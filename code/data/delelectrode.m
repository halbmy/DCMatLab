function [Data1,fim]=delelectrode(Data,fi,removeel)

% DELELECTRODE - Delete measurements with specific electrode
% N = delelectrode(N,fi)

if nargin<2, return; end
if nargin<3, removeel=1; end
if islogical(fi), fi=find(fi); end
la=find(ismember(Data.a,fi));
lb=find(ismember(Data.b,fi));
lm=find(ismember(Data.m,fi));
ln=find(ismember(Data.n,fi));
fim=unique([la(:);lb(:);lm(:);ln(:)]);
Data1=delmeasurement(Data,fim);
if removeel, Data1=deldeadelecs(Data1); end