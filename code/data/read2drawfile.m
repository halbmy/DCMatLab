function N = read2drawfile(datafile,fmt)
% READ2DRAWFILE - read 2d raw file
% N = read2drawfile(filename,formatstring)
% formatstring = [ headerlines n_a n_b n_m n_n rhoa ]

if nargin<2, fmt=[1 1 3 5 7 12]; end
DATA=textread(datafile,'','headerlines',fmt(1));
data=size(DATA,1);
AA=DATA(:,fmt(2));
if fmt(3)>0, BB=DATA(:,fmt(3)); else 
    BB=[]; end
MM=DATA(:,fmt(4));
if fmt(5)>0, NN=DATA(:,fmt(5)); else
    NN=zeros(size(MM)); end
% N.r=DATA(:,max(find(DATA(1,:))));
N.r=DATA(:,fmt(6));
fak=1e3;
AA=round(AA*fak)/fak;BB=round(BB*fak)/fak;
MM=round(MM*fak)/fak;NN=round(NN*fak)/fak;
% message(strcat(mess1,' --> ',titel));
% message(sprintf('%s min=%.1f,max=%.1f',mess2,min(N.r),max(N.r)));
N.elec=unique(sortrows([AA;BB;MM;NN]),'rows');
if size(N.elec,2)<2, N.elec(:,2)=0; end
N.a=zeros(data,1);
N.b=zeros(data,1);
N.m=zeros(data,1);
N.n=zeros(data,1);
N.k=zeros(data,1);

fac=111111;
qel=N.elec(:,1)+N.elec(:,2)*fac;
qa=AA(:,1);qm=MM(:,1);
if size(AA,2)>1, qa=qa+AA(:,2)*fac; end
if size(MM,2)>1, qm=qm+MM(:,2)*fac; end
[aa,bb]=meshgrid(qa,qel);
[ii,jj]=find((aa-bb)==0);
N.a(jj)=ii;
if ~isempty(BB),
    qb=BB(:,1);
    if size(BB,2)>1, qb=qb+BB(:,2)*fac; end
    [aa,bb]=meshgrid(qb,qel);
    [ii,jj]=find((aa-bb)==0);
    N.b(jj)=ii;
end
[aa,bb]=meshgrid(qm,qel);
[ii,jj]=find((aa-bb)==0);
N.m(jj)=ii;
if ~isempty(NN),
    qn=NN(:,1);
    if size(NN,2)>1, qn=qn+NN(:,2)*fac; end
    [aa,bb]=meshgrid(qn,qel);
    [ii,jj]=find((aa-bb)==0);
    N.n(jj)=ii;
end

N.k=getkonf(N);
