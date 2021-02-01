function N=getpseudos(N,dirs)

% GETPSEUDOS - Extract profiles from 3D data set
% N = getpseudos(N);
% adds the following cells to the data structure N
% zweid - the profile data structures
% names - the names of the profiles (e.g. x=-2)
% nr    - the corresponding numbers in the 3D data
% getpseudos(N,1/2) gets only profiles in x/y direction
% getpseudos(N,minnr) gets only profiles with at least minnr electrodes

minnr=1;
if nargin<2, 
    dirs=1:2; % both x and y
else
    if dirs(1)>2, minnr=dirs(1);dirs=1:2; end % minimum number
end
bm=(1:size(N.elec,1))'; % backmap vector
xy='xy';
pro=0;
N.zweid={};N.nr={};N.names={};N.points={};
for i=dirs,  % x/y
    xx=unique(N.elec(:,i)); % all x/y-lines
    for nx=1:length(xx), % all profiles
        aa=(N.elec(N.a,i)==xx(nx))&(N.elec(N.m,i)==xx(nx)); %% Ax=Mx=x
%         ab=ones(size(aa))*(-999);mn=ab+1;
        fnb=find((N.b>0)&(N.n>0));
        fb=find(N.b);fn=find(N.n);
        if fb, aa(fb)=aa(fb)&((N.elec(N.b(fb),i)==xx(nx))); end % Bx=x|B=0
        if fn, aa(fn)=aa(fn)&((N.elec(N.n(fn),i)==xx(nx))); end % Nx=x|N=0
        if fnb, aa(fnb)=aa(fnb)|(N.elec(N.a,i)==xx(nx)&(N.elec(N.b(fnb),i)==N.elec(N.n(fnb),i))); end %Bx=Nx|B=N=0
%         aa(fi)=aa(fi)&((N.elec(N.b(fi),i)==xx(nx))|(N.elec(N.b(fi),i)==N.elec(N.n(fi),i)));
%         aa(fi)=aa(fi)&((N.elec(N.n(fi),i)==xx(nx))|(N.elec(N.b(fi),i)==N.elec(N.n(fi),i)));
%         ab(fi)=N.elec(N.a(fi),i)+N.elec(N.b(fi),i);
%         mn(fi)=N.elec(N.m(fi),i)+N.elec(N.n(fi),i);
        fx=find(aa);%|(ab==mn));
        if ~isempty(fx), % found some on profile
            fel=find(N.elec(:,i)==xx(nx));
            if length(fel)<minnr, break; end
            pro=pro+1;
%             nn.elec=massbandkorr(N.elec(fel,[3-i 3]));
% dont't know why I introduced the next line
%             nn.elec=[0;cumsum(sqrt(sum(diff(N.elec(fel,3-i:3)).^2,2)))]+N.elec(fel(1),3-i);
            nn.elec=N.elec(fel,3-i);
            di=round(diff(sort(nn.elec))*20)/20;
            del=round(median(di)*20)/20;
            nnel=length(nn.elec);
            if length(unique(di))>1, 
                nn.elec=(0:nnel-1)'*del+ceil(nn.elec(1,1)/del)*del; 
            else
                nn.elec=round(nn.elec/del)*del;%ceil
            end
            nn.elec(:,2)=0;%!!!
            bm(:)=0;
            bm(fel)=(1:length(fel))';
            nn.a=bm(N.a(fx(:)));
            nn.m=bm(N.m(fx(:)));
            nn.b=zeros(size(nn.a));
            nn.n=zeros(size(nn.a));
            du=N.b(fx);fi=find(du);
            nn.b(fi)=bm(du(fi));
            du=N.n(fx);fi=find(du);
            nn.n(fi)=bm(du(fi));
            nn.b(nn.b>nnel)=0;
            nn.n(nn.n>nnel)=0;
            el=ones(size(fx))*eps;
            fi=find(N.b(fx));
            el(fi)=N.elec(N.b(fx(fi)),i);
            nn.b(el~=xx(nx))=0;
            fi=find(N.n(fx));
            el(fi)=N.elec(N.n(fx(fi)),i);
            nn.n(el~=xx(nx))=0;            
%             nn.b(N.elec(N.b(fx),i)~=xx(nx))=0; % buggy for N.b==0
%             nn.n(N.elec(N.n(fx),i)~=xx(nx))=0;
            if isfield(N,'r'), nn.r=N.r(fx); end
            if isfield(N,'k'), nn.k=N.k(fx); end
            if isfield(N,'err')&&(length(N.err)>=max(fx)), nn.err=N.err(fx); end
            if isfield(N,'ip')&&(length(N.ip)>=max(fx)), nn.ip=N.ip(fx); end
            if isfield(N,'i')&&(length(N.i)>=max(fx)), nn.i=N.i(fx); end
            if isfield(N,'u')&&(length(N.u)>=max(fx)), nn.u=N.u(fx); end
            name=[xy(i) ' = ' num2str(xx(nx))];
            N.zweid{pro}=nn;
            N.nr{pro}=fx;
            N.names{pro}=name;
            poi(i)=xx(nx);poi(i+2)=xx(nx);
            poi(3-i)=min(nn.elec(:,1));
            poi(5-i)=max(nn.elec(:,1));
            N.points{pro}=poi;
%             showdata2d(nn);
%             text(1,1,name);
%             pause(1.0);
        end
    end
end
if isempty(N.zweid),
    N=rmfield(N,{'zweid','names','nr','points'});
end