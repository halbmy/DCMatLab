function A=showasimage(Data,field,fi,pos,cbarvert)

% SHOWASIMAGE(Data[,field[,indices])

if nargin<5, cbarvert=0; end
if (nargin<2)||isempty(field), 
    if isfield(Data,'t') field=Data.t; end
    if isfield(Data,'r'), field=Data.r; end    
end
if (nargin<3)||isempty(fi), fi=1:length(field); end
if nargin<4, pos=[]; end
if islogical(fi), fi=find(fi); end
if (nargin==3)&&any(fi-round(fi)),
    pos=fi;
    fi=1:length(field);
end
if isstr(Data), Data=readdatafile(Data); end
if min(fi)<1, fi=find(fi); end

nel=1;
if isfield(Data,'elec'), nel=size(Data.elec,1); end
if isfield(Data,'pos'), nel=size(Data.pos,1); end
isalpha=1;
if isalpha, 
    A=ones(nel,nel)*NaN; 
else 
    A=zeros(nel,nel); 
end
if isfield(Data,'a')&&isfield(Data,'m'), %DC case
    for j=1:length(fi), i=fi(j);
        A(Data.a(i),Data.m(i))=field(i); end
    xl='M electrode';yl='A electrode';
else
    xl='';yl='';
end
if isfield(Data,'s')&&isfield(Data,'g'), %TT case
    for j=1:length(fi), i=fi(j);
        A(Data.s(i),Data.g(i))=field(i); end
    yl='shot';xl='geophone';
end
x=1:size(A,2);y=1:size(A,1);
if (nargin<4)&&isfield(Data,'elec'), pos=Data.elec(:,1); end
if (nargin<4)&&isfield(Data,'pos'), pos=Data.pos(:,1); end
if ~isempty(pos), % position for plotting given or auto
    if min(size(pos))>1, % x and z given
        fi0=find(pos(:,end)==0);
        pos(:,end)=-abs(pos(:,end));
        pos(fi0,end)=sum(pos(fi0,1:end-1),2);
        pos(:,1:end-1)=[];
    end
    x=pos(x);y=pos(y); 
    xl=[xl ' pos'];yl=[yl ' pos'];
else
    xl=[xl ' num'];yl=[yl ' num'];
end
fx=find(prod(double(isnan(A))));
fy=find(prod(double(isnan(A)),2));
A(fy,:)=[];A(:,fx)=[];x(fx)=[];y(fy)=[];
imagesc(A);axis ij equal tight;grid on
xt=get(gca,'XTick');yt=get(gca,'YTick');
set(gca,'XTickLabelMode','manual','YTickLabelMode','manual',...
    'XTickLabel',num2strcell(rndig(x(xt),3)),'YTickLabel',num2strcell(rndig(y(yt),3)));

xlabel(xl);ylabel(yl);
if min(field(fi))<0,
   caxis([-1 1]*std(field(fi))*2);
else
   caxis(interperc(field(fi))); 
end
% colormap(b2r);
if cbarvert, colorbar; else colorbar horiz; end
if isalpha, alpha(1-isnan(A)); end
if nargout==0, A=1; end