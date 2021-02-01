function N=createsxhdata(nsel,x0,dx,holes,nbel,rdip,tdip)

% CREATESXHDATA - Create Surface-cross-borehole dataset
% N=createsxhdata(nsel,x0,dx,holes,nbel,rdip,tdip)
% nsel..number of surface electrodes at x0,x0+dx,...
% holes..vector of borehole positions
% nbel..number of borehole electrodes (may be array)
% createsxhdata(25,-12,1,[-6 6],9)
% rdip..receiver is dipole (otherwise pole)
% tdip..transmittor is dipole (otherwise pole)

if nargin<5, error('Too less input arguments'); end
if nargin<6, rdip=0; end
if nargin<7, tdip=0; end

N.elec=(0:nsel-1)'*dx+x0;
N.elec(:,2)=0;
nholes=length(holes);
for n=1:nholes,
    last=nsel+(n-1)*nbel;
    N.elec(last+1:last+nbel,1)=holes(n);
    N.elec(last+1:last+nbel,2)=(1:nbel)'*dx;
end
N.a=[];N.b=[];N.m=[];N.n=[];N.r=[];
% Surface to borehole
for n=1:nholes,
   ib=min(length(nbel),n);
   for s=1:nsel-tdip,
      for b=1:nbel(ib)-rdip,
        N.a(end+1)=s;
        N.m(end+1)=nsel+(n-1)*nbel(ib)+b;
        N.b(end+1)=(N.a(end)+1)*tdip;
        N.n(end+1)=(N.m(end)+1)*rdip;
      end
   end
end
% cross-hole
if nsel==0, % only if no surface electrodes given
    for n=1:nholes,
        ib=min(length(nbel),n);
        for n2=n+1:nholes,
            for b=1:nbel(ib)-tdip,
                for b2=1:nbel(ib)-rdip,
                    N.a(end+1)=nsel+(n-1)*nbel(ib)+b;
                    N.m(end+1)=nsel+(n2-1)*nbel(ib)+b2;
                    N.b(end+1)=(N.a(end)+1)*tdip;
                    N.n(end+1)=(N.m(end)+1)*rdip;
                end
            end
        end
    end
end
N.a=N.a(:);N.b=N.b(:);N.m=N.m(:);N.n=N.n(:);
N.r=ones(size(N.a))*100;
N.k=getkonf(N);