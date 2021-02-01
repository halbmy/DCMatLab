function [fi1,fi2,rez]=getrez(N,field)

% GETREZ - Get reciprocal measurement indices
% [fi1,fi2] = getrez(N);

aa=max(N.a,N.b);
bb=min(N.a,N.b);
mm=max(N.m,N.n);
nn=min(N.m,N.n);
abmn=aa*999999+bb*9999+mm*99+nn;
mnab=mm*999999+nn*9999+aa*99+bb;
fia=find(aa<mm);
fim=find(aa>mm);
[CC,i1,i2]=intersect(abmn(fia),mnab(fim));
fi1=fia(i1);
fi2=fim(i2);
% [N.a(fi1(1:10:100+100)) N.m(fi2(1:10:100+100))]
if nargout~=2,
    if nargin<2, field=N.r; end
    R1=field(fi1);
    R2=field(fi2);
    rezi=(R1-R2)./(R1+R2);
    rez=N.a*0;
    rez(fi1)=rezi;
    rez(fi2)=rezi;
    if nargout<1,
        figure(1);
        plot(R1,R2,'.');xlabel('normal');ylabel('reverse');
        figure(2);
        hist(rez*100,100);xlabel('reciprocity/%');ylabel('frequency');
    end
    std_max_rez=[std(rez) max(abs(rez))]*100; % in %
end