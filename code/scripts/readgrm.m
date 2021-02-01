% load layer
fname='kon2.grm';
fid=fopen(fname);
% C=textscan(fid,'%s','commentStyle','#','delimiter','EOL');c=C{1};
zeile=fgetl(fid);
shotnr=0;
while ~isequal(zeile,-1),
    if findstr(zeile,'SHOT'), zeile=fgetl(zeile); end
    zeile=fgetl(fid);
end
fclose(fid);
