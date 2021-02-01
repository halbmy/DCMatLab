function showrays(Mesh,Shot,ishot,vel,posround)

% SHOWRAYS - Show rays from individual shot point
% showrays(Mesh,Shot,ishot)

if nargin<2, error('Mesh and Shot must be provided!'); end
if nargin<3, ishot=floor(rand*length(Shot.loc))+1; end
if nargin>3, Mesh.cellattr=1./vel; end
if nargin<5, posround=100; end
Mesh.edgelength=sqrt(sum((Mesh.node(Mesh.bound(:,1),:)-Mesh.node(Mesh.bound(:,2),:)).^2,2));
Mesh.edgeslowness=ones(size(Mesh.edgelength));
for i=1:Mesh.nbounds,
    lr=[Mesh.boundleft(i) Mesh.boundright(i)];
    lr(lr==0)=[];
    Mesh.edgeslowness(i)=min(Mesh.cellattr(lr));
end
Mesh.edgetime=Mesh.edgelength.*Mesh.edgeslowness;
l=0;
bmi=min(Mesh.bound,[],2);
bma=max(Mesh.bound,[],2);
aa=round(Shot.pos*posround)/posround;
bb=round(Mesh.node*posround)/posround;
if size(aa,2)<2, aa(:,2)=0; end
% [cc,is]=intersect(bb,aa,'Rows');
[cc,is]=ismember(aa,bb,'Rows');
i=ishot;
start=is(Shot.ns{i});
[dist,prec]=dijkstra(Mesh,start);
for j=1:length(Shot.nx{i}),
    l=l+1;
    w=is(Shot.nx{i}(j));
    while w(end)~=start,
        w=[w prec(w(end))];
        fb=find((bmi==min(w(end-1:end)))&(bma==max(w(end-1:end))));
        cc=[Mesh.boundleft(fb) Mesh.boundright(fb)];cc=cc(1,:);
        cc(cc==0)=[];
        sls=Mesh.cellattr(cc);cc(sls>mean(sls))=[];
    end
    w=fliplr(w);
%     hold on;plot(Mesh.node(w,1),Mesh.node(w,2),'k-','LineWidth',1);hold off;
    set(line(Mesh.node(w,1),Mesh.node(w,2)),'Color','black');
end
