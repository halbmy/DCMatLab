function Shot = pick2shot(filenr,shotpos,geopos)

% PICK2SHOT - Reads (su) pick files into Shot
% Shot = pick2shot(filenr,shotpos,geopos)
% filenr...pickfile without .pick
% shotpos..shot position vector
% geopos...geophone position vector
% Example with 9 shots and 48 electrodes
% filenr=1001:1009; % 9 shots named 1001.pick,...,1009.pick
% shotpos=-1:12:95; % every 12m starting from -1
% geopos=0:2:94; % equdistant

Shot=[];
if length(filenr)~=length(shotpos), 
    error('Size mismatch!');return; 
end
ngeo=length(geopos);
Shot.t=[];
Shot.pos=[geopos(:);shotpos(:)];Shot.pos(:,2)=0;
for i=1:length(filenr),
   Shot.ns{i}=i+ngeo;
   A=load([num2str(filenr(i)) '.pick']);
   Shot.tt{i}=A(:,1)*1000;
   Shot.t=[Shot.t;A(:,1)];
   Shot.nx{i}=A(:,2);
end
%%
[Shot.pos,II]=sortrows(Shot.pos);
AI=zeros(size(II));AI(II)=1:length(II);
for i=1:length(Shot.ns),
   Shot.ns{i}=AI(Shot.ns{i});
   Shot.nx{i}=AI(Shot.nx{i});
end
if nargout<1, plotshot(Shot); end
% savesgtfile('sportplatz.sgt',Shot);