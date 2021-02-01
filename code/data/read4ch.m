function N=read4ch(datfile)
% datfile='ddf3d4.org';
[a,b,m1,n1,m2,n2,m3,n3,m4,n4]=textread(datfile,'%d %d %d %d %d %d %d %d %d %d ','headerlines',2);
n=length(a)*4;
N.a=zeros(n,1);N.b=N.a;N.m=N.a;N.n=N.a;
N.r=ones(n,1);
k=0;
mm=[m1 m2 m3 m4];
nn=[n1 n2 n3 n4];
for i=1:length(a),
    for j=1:4,
        k=k+1;
        N.a(k)=a(i);
        N.b(k)=b(i);
        N.m(k)=mm(i,j);
        N.n(k)=nn(i,j);
    end
end
fi=find(N.m==0);
N.a(fi)=[];N.b(fi)=[];N.m(fi)=[];N.n(fi)=[];N.r(fi)=[];
ma=max([N.a;N.b;N.m;N.n]);
N.elec=(1:ma)'-1;
N.elec(:,2)=0;
message(sprintf('Reading %d data points from 4-channel file %s',length(N.a),datfile));
message(sprintf('Effectivity = %d%%',round(100*length(N.a)/(length(fi)+length(N.a)))));
N.k=getkonf(N);
