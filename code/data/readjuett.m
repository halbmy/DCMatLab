function Shot1 = readjuett(filename,dx)
% READJUETT - Read Juelich type travel time data
% Shot = readjuett(filename)

% filename='D:\Guenther.T\2d\juelich\05_10_12\ht100_13t_cuts.dat';
% filename='D:\Guenther.T\2d\juelich\05_10_12\8-24_13t_all.dat';

dx=5e-3;
if nargin<2, dx=0.0001; end
A=textread(filename,'');
Shot=[];
% Shot.t=A(:,1)*1e-9;
tr=A(:,3:5);
re=A(:,6:8);
tr=round(tr/dx)*dx;re=round(re/dx)*dx;
Shot.pos=unique([tr;re],'rows');fprintf('%d positions found',size(Shot.pos,1));
wi=atan2(Shot.pos(:,2),Shot.pos(:,1));
[swi,ind]=sort(wi);
Shot.pos=Shot.pos(ind,:);
untr=unique(tr,'rows');
[tf,nre]=ismember(re,Shot.pos,'rows');
Shot.t=[];
for i=1:size(untr,1),
    [tf,Shot.ns{i}]=ismember(untr(i,:),Shot.pos,'rows');
    fi=find(ismember(tr,untr(i,:),'rows'));
    Shot.tt{i}=A(fi,1);
    Shot.t=[Shot.t;Shot.tt{i}*1e-9];
    Shot.nx{i}=nre(fi);
end
for i=3:-1:1, if length(unique(Shot.pos(:,i)))==1, Shot.pos(:,i)=[]; end; end
Shot.sd=3e-10;
if size(A,2)>8, %va present
    Shot.va=A(:,9)./Shot.t;
else
    Shot.dist=zeros(size(Shot.t));
    l=0;
    for i=1:length(Shot.ns),
        di=sqrt((Shot.pos(Shot.nx{i},1)-Shot.pos(Shot.ns{i},1)).^2+...
            (Shot.pos(Shot.nx{i},2)-Shot.pos(Shot.ns{i},2)).^2);
        le=length(di);
        Shot.dist(l+1:l+le)=di;
        l=l+le;
    end
    Shot.va=Shot.dist./Shot.t;
end
Shot1=reduceshot(Shot);
fprintf('...reduced to %d\n',size(Shot1.pos,1));
return
clf;plot(Shot.pos(:,1),Shot.pos(:,2),'x');axis equal;grid on;
hold on;ns=cell2mat(Shot.ns);plot(Shot.pos(ns,1),Shot.pos(ns,2),'ro');hold off
return
for i=1:length(Shot.ns), 
    clf;plot(Shot.pos(:,1),Shot.pos(:,2),'.');axis equal;grid on;
    for j=1:length(Shot.nx{i}),
        aa=[Shot.ns{i} Shot.nx{i}(j)];line(Shot.pos(aa,1),Shot.pos(aa,2));
    end; 
    pause
end
