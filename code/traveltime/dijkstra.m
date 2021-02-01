function [dist,prec]=dijkstra(Mesh,start)

% DIJKSTRA - perform dijskstra routing algorithm on mesh
% [dist,prec] = dijkstra(Mesh,start)

U=1:Mesh.nnodes;
dist=ones(Mesh.nnodes,1)*999;dist(start)=0;
prec=zeros(Mesh.nnodes,1);prec(start)=start;
while ~isempty(U),
    [umin,iu]=min(dist(U));
    u=U(iu);U(iu)=[];
    fi1=find(Mesh.bound(:,1)==u);
    fi2=find(Mesh.bound(:,2)==u);fi=[fi1;fi2];
    allv=[Mesh.bound(fi1,2);Mesh.bound(fi2,1)];
%     [fi,fi2]=find(Mesh.bound==u);
%     allv=Mesh.bound(fi,:);allv(allv==u)=[];
    for iv=1:length(allv),
        v=allv(iv);
        newt=dist(u)+Mesh.edgetime(fi(iv));
        if newt<dist(v),
            dist(v)=newt;
            prec(v)=u;
        end
    end
end
% way=[stop prec(stop)];
% while way(end)~=start,
%     way=[way prec(way(end))];
% end
% way=fliplr(way)
% time=dist(stop);
