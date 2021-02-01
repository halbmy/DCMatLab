function savesgtfile(filename,Shot,s,g)

% SAVESGTFILE - Save file to unified data format with Shot,Geophone,Time
% savesgtfile(filename,Shot)

% filename='test.dat';

if nargin<4,
    s=[];g=[];
    for i=1:length(Shot.nx),
        s=[s;ones(length(Shot.nx{i}),1)*Shot.ns{i}];
        g=[g;Shot.nx{i}(:)];
    end
end
nl='\r\n';
sa='%.6f';sa='%g';ss=sa;
fo='#x';ko={'x','y','z'};
for i=2:size(Shot.pos,2), 
    ss=[ss '\t' sa]; 
    fo=[fo '\t' ko{i}];
end
if length(g)==length(Shot.t),
fid=fopen(filename,'w');
fprintf(fid,['%d # shot/geophone points' nl],size(Shot.pos,1));
fprintf(fid,[fo nl]);
fprintf(fid,[ss nl],Shot.pos');
fprintf(fid,['%d # measurements' nl],length(Shot.t));
fprintf(fid,['#s\tg\tt' nl]);
fprintf(fid,['%d\t%d\t%g' nl],[s g Shot.t]');
fclose(fid);
else
  fprintf('Size mismatch: %d positions %d times\n',length(g),length(Shot.t));
end