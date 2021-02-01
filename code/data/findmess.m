function fi=findmess(N,el1,el2)
el1=[el1(:);0];
el2=[el2(:);0];
a1=ismember(N.a,el1);
fi=find(N.b);a1(fi)=a1(fi)&ismember(N.b(fi),el1);
a2=ismember(N.a,el2);
fi=find(N.b);a2(fi)=a2(fi)&ismember(N.b(fi),el2);
m1=ismember(N.m,el1);
fi=find(N.n);m1(fi)=m1(fi)&ismember(N.n(fi),el1);
m2=ismember(N.m,el2);
fi=find(N.n);m2(fi)=m2(fi)&ismember(N.n(fi),el2);
fi=find(((a1+a2).*(m1+m2)>0)+((a1+m1).*(a2+m2)>0)==2);