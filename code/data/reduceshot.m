function Shot1=reduceshot(Shot,dx)

% REDUCESHOT - Reduce shot/geophone numbers by assembling
% Shot1=reduceshot(Shot)

if nargin<2, dx=2e-2; end
di=sqrt(sum(diff([Shot.pos;Shot.pos(1,:)]).^2,2));
npos=size(Shot.pos,1);
fi=find(di<dx);
list=1:npos;
dellist=[];
for i=1:length(fi),
   eins=fi(i);
   zwei=eins+1;
   if zwei>npos, zwei=1; end
   Shot.pos(eins,:)=mean(Shot.pos([eins zwei],:));
   list(zwei)=list(eins);
   dellist(end+1)=zwei;
end
indlist=[0 cumsum(sign(diff(list)))]+1;
Shot1=Shot;
Shot1.pos(dellist,:)=[];
for i=1:length(Shot.ns),
    Shot1.ns{i}=indlist(Shot1.ns{i});
    Shot1.nx{i}=indlist(Shot1.nx{i});
end
% shotimage(Shot1);
di1=sqrt(sum(diff([Shot1.pos;Shot1.pos(1,:)]).^2,2));
fprintf('reduced to %d\n',size(Shot1.pos,1));
if nargout==0,
    plot(Shot.pos(:,1),Shot.pos(:,2),'bo-',Shot1.pos(:,1),Shot1.pos(:,2),'rx-');
end
