% SS=readunifile('sinkhole1_noise1ms.sgt');
SS=readunifile('sinkhole2.sgt');
A=[(1:length(SS.x))' SS.x SS.y SS.z ones(size(SS.x))];
% fid=fopen('sinkhole1.geo','w');
fid=fopen('sinkhole2.geo','w');
fprintf(fid,'%d %.2f %.2f %.2f %d\n',A');
fclose(fid);
B=[(1:length(SS.a))' SS.a SS.m ones(size(SS.a)) SS.t*1000];
% fid=fopen('sinkhole1.dat','w');
fid=fopen('sinkhole2.dat','w');
fprintf(fid,'%d %d %d %d %.2f\n',B');
fclose(fid);