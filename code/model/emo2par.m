function [res,thk] = emo2par(filename)

% EMO2PAR - reads parameters (resistivity&thickness) from emo file
% [res,thk] = emo2par(filename)

% filename='test.emo';
fid=fopen(filename,'r');
zeile=fgetl(fid);
s_ite=sum(('N iterations (NI, 0 -> analysis and forward resp.)'-32).^2);
s_par=sum(('Parameters (0..NIte)'-32).^2);
s_nlay=sum(('Number of layers per model (1..NM)'-32).^2);
iter=0;numlay=0;imodel=0;res={};thk={};
while ischar(zeile),
    if sum((zeile-32).^2)==s_ite, iter=fscanf(fid,'%d'); end
    if sum((zeile-32).^2)==s_nlay, numlay=fscanf(fid,'%d'); end
    zeile=fgetl(fid);
    if sum((zeile-32).^2)==s_par,
        imodel=imodel+1;
        for i=1:iter+1, zeile=fgetl(fid); end
        par=str2num(zeile);
        while par(1)<iter, zeile=fgetl(fid);par=str2num(zeile); end        
        res{imodel}=par(2:numlay+1);
        thk{imodel}=par(numlay+2:2*numlay);
    end
end
if length(res)==1,% only 1 sounding->vector instead of cell
    res=res{1};thk=thk{1};
end
fclose(fid);