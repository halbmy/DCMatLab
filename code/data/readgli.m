function Shot=readgli(filename)

% READGLI - Read GLI arrival time file
% Shot = readgli(filename)
% Columns are:
% SHOT/TRACE marker x   xref    z   t/ms

fid=fopen(filename,'r');
if fid<0, error(['Filename ' filename ' could not be opened!']); end
ishot=0;Shot=[];itrace=0;
zeile=fgetl(fid);
while ischar(zeile),
    typ=sscanf(zeile,'%s5%f');
    nu=sscanf(zeile,'%*s%f%f%f%f%f');
    if length(nu)<5, break; end
    if isequal(typ,'SHOT'),
        if ishot>0, fprintf('%d traces\n',itrace); end
        ishot=ishot+1;
        itrace=0;
        Shot.loc(ishot)=nu(2)-nu(3);
        fprintf('Shot %d at x=%g ',ishot,Shot.loc(ishot));
        Shot.locz(ishot)=nu(4);
        Shot.x{ishot}=[];
        Shot.z{ishot}=[];
        Shot.tt{ishot}=[];
    elseif isequal(typ,'TRACE'),
        itrace=itrace+1;
        Shot.x{ishot}(end+1)=nu(2)-nu(3);
        Shot.z{ishot}(end+1)=nu(4);
        Shot.tt{ishot}(end+1)=nu(5);
    else
        nu
    end
    zeile=fgetl(fid);
end
fprintf('%d traces\n',itrace);
fclose(fid);
pos=[];Shot.t=[];posxz=[];
for i=1:length(Shot.loc),
    pos=[pos;Shot.x{i}(:);Shot.loc(i)];
    posxz=[posxz;Shot.x{i}(:) Shot.z{i}(:)];
    Shot.t=[Shot.t;Shot.tt{i}(:)];
end
% pos=round(pos*1000)/1000;posxz=round(posxz*1000)/1000;
Shot.pos=unique(pos,'rows');
posxz=unique(posxz,'rows');
posxz(find(diff(posxz(:,1))==0),:)=[];
Shot.pos(:,2)=round(interp1(posxz(:,1),posxz(:,2),Shot.pos(:,1),'linear','extrap')*100)/100;
for i=1:length(Shot.loc),
    Shot.ns{i}=find(Shot.pos(:,1)==Shot.loc(i));
    [C,ia]=intersect(Shot.pos(:,1),Shot.x{i});
    Shot.nx{i}=ia;
end