% [fname,fpath]=uigetfile('*.gem');
gemfile=fullfile(fpath,fname);
N=readgemfile(gemfile);
fi1=find((N.r>3000)|(N.r<100)|(abs(N.k)>10000));
N.a(fi1)=[];N.b(fi1)=[];N.m(fi1)=[];N.n(fi1)=[];N.k(fi1)=[];N.r(fi1)=[];
for ii=1:3,
    am=min(N.a,N.m)*1000+max(N.a,N.m); % i.A. 4-Punkt (AB + MN sortieren, dann AM sortieren, dann a*100000+...)
    [uam,i,j]=unique(am);[sj,ind]=sort(j);di=diff(sj,2);
    fi=find(di==0);sj(fi)=[];ind(fi)=[];
    rez=2*abs(N.r(ind(1:2:end-1))-N.r(ind(2:2:end)))./(N.r(ind(1:2:end-1))+N.r(ind(2:2:end)));
    fi2=find(rez>0.1);
    fi1=[ind(fi2*2);ind(fi2*2-1)];
    N.a(fi1)=[];N.b(fi1)=[];N.m(fi1)=[];N.n(fi1)=[];N.k(fi1)=[];N.r(fi1)=[];
end
plot(N.r(ind(1:2:end-1)),N.r(ind(2:2:end)),'.');
