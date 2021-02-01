function N=combdata3d(N1,N2)

% COMBDATA3D - Combine 3d data sets
% N=combdata3d(N1,N2)

l1=length(N1.elec);n1=length(N1.a);
l2=length(N2.elec);n2=length(N2.a);
[N.elec,I,J]=unique([N1.elec;N2.elec],'rows');
N.a=zeros(n1+n2,1);
N.b=N.a;N.m=N.a;N.n=N.a;
if isfield(N1,'r')&&isfield(N2,'r'), N.r=[N1.r(:);N2.r(:)]; end
if isfield(N1,'rho')&&isfield(N2,'rho'), N.rho=[N1.rho(:);N2.rho(:)]; end
if isfield(N1,'ip')&&isfield(N2,'ip'), N.ip=[N1.ip(:);N2.ip(:)]; end
if isfield(N1,'err')&&isfield(N2,'err'), N.err=[N1.err(:);N2.err(:)]; end
if isfield(N1,'i')&&isfield(N2,'i'), N.i=[N1.i(:);N2.i(:)]; end
if isfield(N1,'u')&&isfield(N2,'u'), N.u=[N1.u(:);N2.u(:)]; end
N.k=[N1.k(:);N2.k(:)];
N.a(1:n1)=J(N1.a);
N.a(n1+1:n1+n2)=J(N2.a+l1);
N.m(1:n1)=J(N1.m);
N.m(n1+1:n1+n2)=J(N2.m+l1);
fi=find(N1.b);
N.b(fi)=J(N1.b(fi));
fi=find(N2.b);
N.b(fi+n1)=J(N2.b(fi)+l1);
fi=find(N1.n);
N.n(fi)=J(N1.n(fi));
fi=find(N2.n);
N.n(fi+n1)=J(N2.n(fi)+l1);