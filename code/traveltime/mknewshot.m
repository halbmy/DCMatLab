function Shot=mknewshot(npos,dshot,dx)

% MKNEWSHOT - Make new shot
% Shot = mknewshot(npos,dshot,dx)
%        pos  - number of shot/geophone positions
%               or positions itself
%        dshot - shot every dshot positions
%        dx    - geophone distance

if nargin<1, npos=41; end % default values
if nargin<2, dshot=5; end
if nargin<3, dx=1; end
if length(npos)>1, % positions given
    Shot.pos=npos;
    npos=size(Shot.pos,1);
else % number of 
    Shot.pos=(0:npos-1)';Shot.pos(1,2)=0;
end
Shot.t=[];
ishot=0;nshot=1-dshot; % counters
while(nshot<npos), % shot while not at end
    ishot=ishot+1;nshot=nshot+dshot;
    if nshot>npos, nshot=npos; end
    Shot.ns{ishot}=nshot;
    aa=(1:npos)';aa(nshot)=[];Shot.nx{ishot}=aa;
    Shot.loc(ishot)=Shot.pos(nshot,1); %historical
    Shot.x{ishot}=Shot.pos(aa); % only for plotting
    Shot.tt{ishot}=abs(aa-nshot);
    Shot.t=[Shot.t;Shot.tt{ishot}/1000];
end
if nargout<1, plotshot(Shot); end