function Shot = readtom(filename)

% READTOM - Read TOM travel time file
% Shot = readtom(filename)

Shot=[];
% SGT=sortrows(load(filename),[3 6]);%x of shot and receiver
fid=fopen(filename,'r');
SGT=fscanf(fid,'%f',[8 Inf])';
fclose(fid);
Shot.pos=unique([SGT(:,[3 5]);SGT(:,[6 8])],'rows');
Shot.t=SGT(:,1);
ishot=0;ges=0;
while ~isempty(SGT),
    la=min(find(SGT(:,3)~=SGT(1,3)))-1;
    dt=0;
    if isempty(la), la=size(SGT,1); end
    ishot=ishot+1;
    Shot.loc(ishot)=SGT(1,3);Shot.locz(ishot)=SGT(1,5);
    Shot.x{ishot}=SGT(1:la,6);Shot.z{ishot}=SGT(1:la,8);
    self=find(Shot.x{ishot}==Shot.loc(ishot));
    if ~isempty(self),
        dt=min(SGT(self,1));
        Shot.x{ishot}(self)=[];Shot.z{ishot}(self)=[];
        Shot.t(ges+self)=[];SGT(self,:)=[];la=la-length(self);
        Shot.t(ges+1:ges+la)=Shot.t(ges+1:ges+la)-dt;
    end
    [aa,Shot.ns{ishot}]=ismember(Shot.loc(ishot),Shot.pos(:,1)); %rows
    [aa,Shot.nx{ishot}]=ismember(Shot.x{ishot},Shot.pos(:,1)); %rows
    Shot.tt{ishot}=SGT(1:la,1)-dt;
    SGT(1:la,:)=[];ges=ges+la;
end