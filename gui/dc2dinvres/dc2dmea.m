function MEA = dc2dmea(x,z,sigmac,sx,sz,FOR)
% DC2DMEA - DC 2D Forward calculation for multielectrode
% MEA=dc2dmea(x,z,sigmac,ex[,ez])
% x,z..Grid nodes, sigmac..cell conductivities
% ex,ez..electrode positions
% MEA..Matrix of Multielectrode potentials

if nargin<4, error('Usage: mea = dc2dmea(x,z,sigmac,sx,sz)'); end
if nargin<5, sz=zeros(size(sx)); end
if nargin<6, FOR=struct('acc',1e-4,'tol',1e-4,'maxit',200); end
if ~isfield(FOR,'acc'), FOR.acc=1e-4; end
if ~isfield(FOR,'tol'), FOR.tol=1e-4; end
if ~isfield(FOR,'maxit'), FOR.maxit=200; end
smax = length(sx);
XMAX = length(x);
ZMAX = length(z);

sq=sigmac(1,1);

% boundary conductivities
sigmac = [zeros(1,XMAX+1)
    [zeros(ZMAX-1,1) sigmac zeros(ZMAX-1,1)]
    zeros(1,XMAX+1)];

% grid sizes
dx = abs(diff(x));
dx = [dx(1); dx; dx(XMAX - 1)];
dz = abs(diff(z));
dz = [dz(1); dz; dz(ZMAX - 1)];
x = repmat(x,ZMAX,1);
z = reshape(repmat(z',XMAX,1),(XMAX)*(ZMAX),1);

% wavenumbers
wmax = 16;
ing = 8;
inl = wmax - ing;
k0 = 1/2/min(dx);

% Secondary potential
% v=zeros(XMAX*ZMAX,1);
% Total potential
phi=zeros(XMAX*ZMAX,smax);
% phis=phi;

% Gauss-Legendre abszissas ks and weights ws, scaling to real wavenumbers 
% 1 ... ing
[ks, ws] = gauleg(0, 1, ing);
wn = ks .* ks * k0;
weight = 2 * ws .* ks * k0 / pi;

% Gauss-Laguerre abszissas ks and weights ws, scaling to real wavenumbers 
% ing+1 ... wmax
[ks, ws] = gaulag(inl);
wn = [wn, k0 * (ks + 1)];
weight = [weight, k0 * exp(ks) .* ws / pi];

% conductivities and distances needed for coupling coefficients
s1 = reshape(sigmac(2:ZMAX+1, 2:XMAX+1)',[XMAX*ZMAX,1]);
s2 = reshape(sigmac(2:ZMAX+1,1:XMAX)',[XMAX*ZMAX,1]);
s3 = reshape(sigmac(1:ZMAX, 1:XMAX)',[XMAX*ZMAX,1]);
s4 = reshape(sigmac(1:ZMAX, 2:XMAX+1)',[XMAX*ZMAX,1]);
d1 = repmat(dx(2:XMAX+1),ZMAX,1);
d2 = reshape(repmat(dz(2:ZMAX+1)',XMAX,1),XMAX*ZMAX,1);
d3 = repmat(dx(1:XMAX),ZMAX,1);
d4 = reshape(repmat(dz(1:ZMAX)',XMAX,1),XMAX*ZMAX,1);

% coupling coefficients for A(sigmac)
% A0 = -1/2*[ce cs cw cn cp+ce+cs+cw+cn]
A0 = -1/2*[(d4.*s4 + d2.*s1)./d1,...
        (d1.*s1 + d3.*s2)./d2,...
        (d2.*s2 + d4.*s3)./d3,...
        (d3.*s3 + d1.*s4)./d4,...
        -(d2.*(d1.*s1 + d3.*s2) + d4.*(d3.*s3 + d1.*s4))/2];

sigma_sq = zeros(ZMAX+1,XMAX+1);

sigma_sq(2:ZMAX,2:XMAX) = sq - sigmac(2:ZMAX,2:XMAX);
s1 = reshape(sigma_sq(2:ZMAX+1, 2:XMAX+1)',[XMAX*ZMAX,1]);
s2 = reshape(sigma_sq(2:ZMAX+1,1:XMAX)',[XMAX*ZMAX,1]);
s3 = reshape(sigma_sq(1:ZMAX, 1:XMAX)',[XMAX*ZMAX,1]);
s4 = reshape(sigma_sq(1:ZMAX, 2:XMAX+1)',[XMAX*ZMAX,1]);
% coupling coefficients for A(sq - sigmac)
% A0_sq = -1/2*[ce cs cw cn cp+ce+cs+cw+cn]
A0_sq = -1/2*[(d4.*s4 + d2.*s1)./d1,...
        (d1.*s1 + d3.*s2)./d2,...
        (d2.*s2 + d4.*s3)./d3,...
        (d3.*s3 + d1.*s4)./d4,...
        -(d2.*(d1.*s1 + d3.*s2) + d4.*(d3.*s3 + d1.*s4))/2];

tic
h = waitbar(0, '2D DC Forward running ...');
for waveindex = 1:wmax
    waitbar(waveindex/wmax,h);
    % assembly of sparse matrix A(sigma)
    r = sqrt((x(fix(XMAX/2)) - x).^2 + (z(1) - z).^2);
    bound = mixed(A0,XMAX,ZMAX,x,z,x(fix(XMAX/2)),z(1),r,wn(waveindex));
    % bound = dirichlet(A0,XMAX,ZMAX);
    ss=sum(A0(:,1:4),2);
    
    AA = spdiags([A0(:,1:4),...
            (-ss + bound + wn(waveindex)^2*A0(:,5))],...
        [-1, -XMAX, 1, XMAX, 0],...
        XMAX*ZMAX,XMAX*ZMAX)';
    
    % AA = spdiags(diag(AA) + bound, 0, AA);

    R=cholinc(AA,FOR.tol);
    
    % assembly of sparse matrix A(sq-sigma)
    A = spdiags([A0_sq(:,1:4),...
            (-sum(A0_sq(:,1:4),2) + wn(waveindex)^2*A0_sq(:,5))],...
        [-1, -XMAX, 1, XMAX, 0],...
        XMAX*ZMAX,XMAX*ZMAX)';
    for is = 1:smax
        % radii sources-gridnodes
        r = sqrt((sx(is) - x).^2 + (sz(is) - z).^2);
        ind0 = find(~r);
        ind1 = find(r);
        
        % analytical potential for wavenumber wn(waveindex)
        v0(ind1) = 1/(pi*sq) * besselmx(real('K'),0,wn(waveindex)*r(ind1),0);
        v0(ind0) = 0;
        v0 = v0(:);
        
        % calculation of right side b
        b = A*v0;       
        % additional boundary terms
        bound = neumann(A0_sq,XMAX,ZMAX,x,z,sx(is),sz(is),r,sq,wn(waveindex));
        b = b + bound;
        
        % additional boundary terms
        % bound = dirichlet(A0,XMAX,ZMAX);
        % b = b + bound;
        [v,flag,err,iter]=pcg(AA,b,FOR.acc,FOR.maxit,R',R);
        %fprintf('Last error %g at %d iterations.\n',err,iter);
        % Integrate over wavenumbers
        phi(:,is)=phi(:,is)+v*weight(waveindex);        
    end
end

% phis=phi;

%Add normal potential in x-z-domain
for is=1:smax,
    r = sqrt((sx(is) - x).^2 + (sz(is) - z).^2);
    ind0 = find(~r);
    ind1 = find(r);    
    v0(ind1) = 1/(2*pi*sq)./r(ind1);
    v0(ind0) = max(v0(ind1));
    phi(:,is)=phi(:,is)+v0(:);
end

runtime=toc;
close(h);

x = x(1:XMAX);
z = z(XMAX:XMAX:ZMAX*XMAX);
phi=permute(reshape(phi,XMAX,ZMAX,smax),[2 1 3]);
%phis=permute(reshape(phis,XMAX,ZMAX,smax),[2 1 3]);
MEA=zeros(smax,smax);

for is = 1:smax,
  MEA(:,is)=interp1(x,squeeze(phi(1,:,is)),sx);
end
message(sprintf('Forward Calculation: runtime=%.1fs =%.1fs/electrode,',runtime,runtime/smax));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction neumann:
% inhomogeneous neumann boundary condition
function bound = neumann(A,XMAX,ZMAX,x,z,sx,sz,r,sigma,lambda)

bound = zeros(XMAX*ZMAX,1);

% left edge
i = (1 : XMAX : (ZMAX-1)*XMAX+1);
d = abs(x(2) - x(1));
bound(i) = bound(i) - A(i,1).*d.*lambda./pi./sigma.*...
    besselmx(real('K'),1,lambda*r(i),0).*abs(sx-x(1))./r(i);

% right edge
i = (XMAX : XMAX : ZMAX*XMAX);
d = abs(x(XMAX) - x(XMAX-1));
bound(i) = bound(i) - A(i,3).*d.*lambda./pi./sigma.*...
    besselmx(real('K'),1,lambda*r(i),0).*abs(sx-x(XMAX))./r(i);

% bottom edge
i = ((ZMAX-1)*XMAX+1 : ZMAX*XMAX);
d = abs(z(XMAX*ZMAX) - x((XMAX-1)*ZMAX));
bound(i) = bound(i) - A(i,4).*d.*lambda./pi./sigma.*...
    besselmx(real('K'),1,lambda*r(i),0).*abs(sz-z(XMAX*ZMAX))./r(i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction dirichlet:
% homogeneous dirichlet boundary condition
function bound = dirichlet(A,XMAX,ZMAX);

bound = zeros(XMAX*ZMAX,1);

% left edge
i = [1 : XMAX : (ZMAX-1)*XMAX+1];
bound(i) = bound(i) - A(i,1);

% right edge
i = [XMAX : XMAX : ZMAX*XMAX];
bound(i) = bound(i) - A(i,3);

% bottom edge
i = [(ZMAX-1)*XMAX+1 : ZMAX*XMAX];
bound(i) = bound(i) - A(i,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunction mixed:
% homogeneous mixed boundary conditions using asymptotical behaviour of
% potential phi ~ K0(lambda*r)
function bound = mixed(A,XMAX,ZMAX,x,z,sx,sz,r,lambda);

bound = zeros(XMAX*ZMAX,1);

% left edge
i = [1 : XMAX : (ZMAX-1)*XMAX+1]';
d = abs(x(2) - x(1));
kr0 = real(besselmx(real('K'),0,lambda*r(i),0));
kr1 = real(besselmx(real('K'),1,lambda*r(i),0));
kr0(find(~kr0)) = 1;
kr1(find(~kr0)) = 1;
bound(i) = bound(i) ...
    - A(i,1)*d*lambda .* (kr1 ./ kr0) .* abs(sx-x(1))./r(i);

% right edge
i = [XMAX : XMAX : ZMAX*XMAX]';
d = abs(x(XMAX) - x(XMAX-1));
kr0 = real(besselmx(real('K'),0,lambda*r(i),0));
kr1 = real(besselmx(real('K'),1,lambda*r(i),0));
kr0(find(~kr0)) = 1;
kr1(find(~kr0)) = 1;
bound(i) = bound(i) ...
    - A(i,3)*d*lambda .* (kr1 ./ kr0) .* abs(sx-x(XMAX))./r(i);

% bottom edge
i = [(ZMAX-1)*XMAX+1 : ZMAX*XMAX]';
d = abs(z(XMAX*ZMAX) - x((XMAX-1)*ZMAX));
kr0 = real(besselmx(real('K'),0,lambda*r(i),0));
kr1 = real(besselmx(real('K'),1,lambda*r(i),0));
kr0(find(~kr0)) = 1;
kr1(find(~kr0)) = 1;
bound(i) = bound(i) ...
    - A(i,4)*d*lambda .* (kr1 ./ kr0) .* abs(sz-z(XMAX*ZMAX))./r(i);
