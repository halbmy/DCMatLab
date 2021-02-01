function Shot = readgrmfile(filename)

% READGRM - Read Gremix *.grm file with arrival times
% Shot = readgrm(filename)

fid=fopen(filename,'r');
zeile=fgetl(fid);
isact=0;ishot=0;Shot=[];
while ischar(zeile),
    if isempty(zeile), isact=0; end
    if isact==1,
        aa=sscanf(zeile,'%f');
        if length(aa)>1,
            Shot.x{ishot}(end+1)=aa(1);
            Shot.z{ishot}(end+1)=aa(3);
            Shot.tt{ishot}(end+1)=aa(2);
        else
            isact=0;
        end
    else
        if strfind(zeile,'SHOT'),
            loc=sscanf(zeile,'%*s%*s%*s%f%*s');
            locz=sscanf(zeile,'%*s%*s%*s%*s%*s%f%*s');
            zeile=fgetl(fid);
            isact=1;
            ishot=ishot+1;
            fprintf('Shot %d at x=%.1f\n',ishot,loc);
            Shot.x{ishot}=[];
            Shot.z{ishot}=[];
            Shot.tt{ishot}=[];
            Shot.loc(ishot)=loc;%47-loc;
            Shot.locz(ishot)=locz;
        end
    end
    zeile=fgetl(fid);
end
fclose(fid);
pos=[];Shot.t=[];posxz=[];
for i=1:length(Shot.loc),
    pos=[pos;Shot.x{i}(:);Shot.loc(i)];
    posxz=[posxz;Shot.x{i}(:) Shot.z{i}(:)];
    Shot.t=[Shot.t;Shot.tt{i}(:)];
end
Shot.pos=unique(pos,'rows');
posxz=unique(posxz,'rows');
Shot.pos(:,2)=round(interp1(posxz(:,1),posxz(:,2),Shot.pos(:,1),'linear','extrap')*100)/100;
for i=1:length(Shot.loc),
    Shot.ns{i}=find(Shot.pos(:,1)==Shot.loc(i));
    [C,ia]=intersect(Shot.pos(:,1),Shot.x{i});
    Shot.nx{i}=ia;
end