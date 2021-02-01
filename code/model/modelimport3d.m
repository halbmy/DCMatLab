function [M,x,y,z,Cov]=modelimport3d(filename)

% MODELIMPORT3D - Import grid model from ASCII file
% [M,x,y,z] = modelimport3d(filename)
% Model = modelimport3d(filename)

% fid=fopen(filename,'r'); 
% zeile='';
% while length(zeile)==0, zeile=destrip(fgetl(fid)); end
% XYZR=fscanf(fid,'%f',[7,Inf])';
% XYZR(end+1,:)=sscanf(zeile,'%f',[7,1])';
% fclose(fid);
XYZR=textread(filename,'','commentstyle','shell','headerlines',1);
x=unique(XYZR(:,1));x(end+1)=max(XYZR(:,2));
y=unique(XYZR(:,3));y(end+1)=max(XYZR(:,4));
z=unique(XYZR(:,5));z(end+1)=max(XYZR(:,6));
if size(XYZR,2)<8, XYZR(:,8)=0; end
M=ones(length(x)-1,length(y)-1,length(z)-1);
Cov=zeros(size(M));
[tf,ii]=ismember(XYZR(:,1),x);
[tf,jj]=ismember(XYZR(:,3),y);
[tf,kk]=ismember(XYZR(:,5),z);
for l=1:length(ii),
    M(ii(l),jj(l),kk(l))=XYZR(l,7);
    Cov(ii(l),jj(l),kk(l))=XYZR(l,8);
end
if nargout==1,
    du.M=M;
    M=du;
    M.x=x;
    M.y=y;
    M.z=z;
    M.Bg=0;
end

function zeile=destrip(zeile)
% strip string from comments (with # character)
aa=strfind(zeile,'#');
if ~isempty(aa), zeile=zeile(1:aa(1)-1); end
