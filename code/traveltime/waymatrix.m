function [W,t] = waymatrix(Mesh,Shot,vel,posround,offset)

% WAYMATRIX - Calculate waymatrix
% W = waymatrix(Mesh,Shot,velocity)
% [W,t] = waymatrix(Mesh,Shot,velocity)

if nargin<2, error('Mesh and Shot must be specified!'); end
if (nargin<4)||isempty(posround), posround=1000; end
if nargin>2, Mesh.cellattr=1./vel; end
if isfield(Shot,'ns'),
    nShot=length(Shot.ns);
else
    uns=unique(Shot.s);
    nShot=length(uns);
end
nModel=Mesh.ncells;
Mesh.edgelength=sqrt(sum((Mesh.node(Mesh.bound(:,1),:)-Mesh.node(Mesh.bound(:,2),:)).^2,2));
Mesh.edgeslowness=ones(size(Mesh.edgelength));
for i=1:Mesh.nbounds,
    lr=[Mesh.boundleft(i) Mesh.boundright(i)];
    lr(lr==0)=[];
    Mesh.edgeslowness(i)=min(Mesh.cellattr(lr));
end
Mesh.edgetime=Mesh.edgelength.*Mesh.edgeslowness;
Wsize2=Mesh.ncells;
if (nargin>4)&&(length(offset)==nShot),
    Wsize2=Wsize2+nShot;
end
W=spalloc(length(Shot.t),Wsize2,round(length(Shot.t)*Mesh.ncells/10));
l=0;
bmi=min(Mesh.bound,[],2);
bma=max(Mesh.bound,[],2);
aa=round(Shot.pos*posround)/posround;
bb=round(Mesh.node*posround)/posround;
if size(aa,2)<2, aa(:,2)=0; end
% [cc,is]=intersect(bb,aa,'Rows');
[cc,is]=ismember(aa,bb,'Rows');
if find(cc==0), error(['Shot point not found: ' num2str(find(cc==0)')]); end
if isfield(Shot,'ns'),
    for i=1:length(Shot.ns),
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
                W(l,cc)=Mesh.edgelength(fb(1))/length(cc);
            end
            w=fliplr(w);
            if nargout<1,
                hold on;plot(Mesh.node(w,1),Mesh.node(w,2),'k-','LineWidth',1);hold off;
            end
        end
        if nargout<1, break; end
    end
else
    for i=1:length(uns),
        fi=find(uns==uns(i));
        start=is(uns(i));
        [dist,prec]=dijkstra(Mesh,start);
        for j=1:length(fi),
            l=l+1;
            w=is(Shot.g(fi(j)));
            while w(end)~=start,
                w=[w prec(w(end))];
                fb=find((bmi==min(w(end-1:end)))&(bma==max(w(end-1:end))));
                cc=[Mesh.boundleft(fb) Mesh.boundright(fb)];cc=cc(1,:);
                cc(cc==0)=[];
                sls=Mesh.cellattr(cc);cc(sls>mean(sls))=[];
                W(l,cc)=Mesh.edgelength(fb(1))/length(cc);
            end
            w=fliplr(w);
            if nargout<1,
                hold on;plot(Mesh.node(w,1),Mesh.node(w,2),'k-','LineWidth',1);hold off;
            end
        end
    end
end
if (nargin>4)&&(length(offset)==nShot),
   j=0;
%    nBound=size(C,1);
   for i=1:nShot,
       nRec=length(Shot.nx{i});
       W(j+1:j+nRec,nModel+i)=1;
       j=j+nRec;
%        C(nBound+i,nModel+i)=olam;
   end
   deltaModel(nModel+1:nModel+nShot)=0;
   rhs=[Mesh.cellattr;offset];
   t=W*rhs;
else
    t=W*Mesh.cellattr(:);
end