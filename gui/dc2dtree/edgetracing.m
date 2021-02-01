function [time,way]=edgetracing(way,time,stop)

% EDGDETRACING - Trace ray along edges recursively
% [time,waypoints] = edgetracing(start,0,stop)
% calls itself with current time and waypath

global Mesh mintime
if length(way)>20, time=mintime+1;return; end
fi=[find(Mesh.bound(:,1)==way(end));find(Mesh.bound(:,2)==way(end))];
aa=Mesh.bound(fi,:)';aa=aa(:);aa(aa==way(end))=[]; % find neighbors
ism=ismember(aa,way);aa(ism)=[];fi(ism)=[]; % delete passed points
if isempty(aa), time=mintime+1;return; end % dead end
for i=1:length(aa), % all remaining neighbors
    tims(i)=time+Mesh.edgetime(fi(i)); % add time
    ways{i}=[way aa(i)]; % append node
    if (aa(i)==stop), % arrived at the end
        mintime=min(tims(i),mintime); % update guiness book
    else % recursion
        if (tims(i)<mintime), [tims(i),ways{i}]=edgetracing(ways{i},tims(i),stop); else tims(i)=mintime+1; end
    end
end
[time,i]=min(tims); % determine shortest way
way=ways{i}; % return the corresponding raypath