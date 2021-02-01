function [W,t] = waymatrix(Mesh,Shot,vel)

% WAYMATRIX - Calculate waymatrix
% W = waymatrix(Mesh,Shot,velocity)
% [W,t] = waymatrix(Mesh,Shot,velocity)

if nargin<2, error('Mesh and Shot must be specified!'); end
if nargin>2, Mesh.cellattr=1./vel; end
% Mesh.edgelength=sqrt(sum((Mesh.node(Mesh.bound(:,1),:)-Mesh.node(Mesh.bound(:,2),:)).^2,2));
% bmi=min(Mesh.bound,[],2);
% bma=max(Mesh.bound,[],2);
% boundneigh=Mesh.boundneigh;
% boundneigh(boundneigh==0)=Mesh.ncells+1;
% for i=1:Mesh.nbounds,
% %     lr=Mesh.boundneigh;%[Mesh.boundleft(i) Mesh.boundright(i)];
% %     lr(lr==0)=[];
%     Mesh.edgeslowness(i)=min(Mesh.cellattr(boundneigh(i,:)));
% end

% % transferred to loadmesh (dobound)
% Mesh.edge=unique(sort([Mesh.bound(:,1:2);Mesh.bound(:,1:2:3);Mesh.bound(:,2:3)],2),'rows');
% bmi=Mesh.edge(:,1);bma=Mesh.edge(:,2);
% Mesh.edgelength=sqrt(sum((Mesh.node(bmi,:)-Mesh.node(bma,:)).^2,2));
% Mesh.edgeslowness=ones(size(Mesh.edgelength));
% Mesh.edgeneigh=zeros(length(Mesh.edge),5);
% for i=1:length(Mesh.edge),
%     [i1,j1]=find(Mesh.cell==ee(i,1));
%     [i2,j2]=find(Mesh.cell==ee(i,2));
%     fi=intersect(i1,i2);
%     Mesh.edgeneigh(i,1:length(fi))=fi;
% end

edgeneigh=Mesh.edgeneigh;
edgeneigh(edgeneigh==0)=Mesh.ncells+1;
Mesh.cellattr(Mesh.ncells+1)=1e15;
Mesh.edgeslowness=min(Mesh.cellattr(edgeneigh),[],2);
Mesh.edgetime=Mesh.edgelength.*Mesh.edgeslowness;
W=spalloc(length(Shot.t),Mesh.ncells,round(length(Shot.t)*Mesh.ncells/10));
l=0;
aa=round(Shot.pos*1000)/1000;
bb=round(Mesh.node*1000)/1000;
if size(aa,2)<2, aa(:,2)=0; end
% [cc,is]=intersect(bb,aa,'Rows');
[cc,is]=ismember(aa,bb,'Rows');
if find(cc==0), error(['Shot point not found: ' num2str(find(cc==0)')]); end
for i=1:length(Shot.ns),
    start=is(Shot.ns{i});
    [dist,prec]=dijkstra(Mesh,start);
    for j=1:length(Shot.nx{i}),
        l=l+1;
        w=is(Shot.nx{i}(j));
        while w(end)~=start,
            w=[w prec(w(end))];
            fb=find((Mesh.edge(:,1)==min(w(end-1:end)))&(Mesh.edge(:,2)==max(w(end-1:end))));
            cc=edgeneigh(fb,:);%[Mesh.boundleft(fb) Mesh.boundright(fb)];
            cc=cc(1,:);
            cc(cc>Mesh.ncells)=[];
            sls=Mesh.cellattr(cc);cc(sls>mean(sls)+1e-6)=[];
            W(l,cc)=Mesh.edgelength(fb(1))/length(cc);
        end
        w=fliplr(w);
        if nargout<1,
            hold on;plot(Mesh.node(w,1),Mesh.node(w,2),'k-','LineWidth',1);hold off;
        end
    end
    if nargout<1, break; end
end
t=W*Mesh.cellattr(1:Mesh.ncells);