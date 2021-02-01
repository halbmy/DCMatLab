function rec = getshotrec(Shot)

% GETSHOTREC - Get shot reciprocity
% rec = getshotrec(Shot)
% Shot..shot structure with pos,ns,nx and t
% rec..relative reciprocity

l=0;
for i=1:length(Shot.ns),
    for j=1:length(Shot.nx{i}),
        l=l+1;
        nn=sort([Shot.ns{i} Shot.nx{i}(j)]);
        po(l)=nn(1)*10000+nn(2);
    end
end
rec=zeros(l,1);
l=0;
for i=1:length(Shot.ns),
    for j=1:length(Shot.nx{i}),
        l=l+1;
        fi=find(po==po(l));
        if any(fi), 
            tt=Shot.t(fi);
            rec(l)=std(tt)/mean(tt); 
        end
    end
end
