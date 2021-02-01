function W = waymatrix(Mesh)

Mesh.edgelength=sqrt(sum((Mesh.node(Mesh.bound(:,1),:)-Mesh.node(Mesh.bound(:,2),:)).^2,2));
Mesh.edgeslowness=ones(size(Mesh.edgelength));
for i=1:Mesh.nbounds,
    lr=[Mesh.boundleft(i) Mesh.boundright(i)];
    lr(lr==0)=[];
    Mesh.edgeslowness(i)=min(Mesh.cellattr(lr));
end
Mesh.edgetime=Mesh.edgelength.*Mesh.edgeslowness;
fi=find(Mesh.nodemarker==-99);lfi=length(fi);
W=zeros(lfi*(lfi-1)/2,Mesh.ncells);l=0;
bmi=min(Mesh.bound,[],2);
bma=max(Mesh.bound,[],2);
for i=1:lfi,
    start=fi(i);
    [dist,prec]=dijkstra(Mesh,start);
    for j=i+1:lfi,
        l=l+1;
        w=fi(j);
        while w(end)~=start,
            w=[w prec(w(end))];
            fb=find((bmi==min(w(end-1:end)))&(bma==max(w(end-1:end))));
            cc=[Mesh.boundleft(fb) Mesh.boundright(fb)];
            cc(cc==0)=[];
%             sls=Mesh.cellattr(cc);cc(sls>mean(sls))=[];
            W(l,cc)=Mesh.edgelength(fb)/length(cc);
        end
        w=fliplr(w);
    end
end