function [M,x,z,Cov]=modelimport2d(modfile)

% MODELIMPORT2d - Imports 2d model from ASCII-File
% [M,x,z] = modelimport2d(modelfile) or
% Mod = modelimport2d(modelfile)
% Model .. structure with x/z/M fields
% M .. model parameter matrix (length(x)-1)x(length(z)-1)
% x/z .. grid nodes in x/z direction
% the model file is assumed to be written as follows
% x_1 x_2 z_1 z_2 rho # for 1.st model cell
% ...
% x_1 x_2 z_1 z_2 rho #for last model cell

fid=fopen(modfile,'r');
if fid<0, error(['File ' modfile ' does not exist!']); end
zeile='';
while isempty(zeile), zeile=destrip(fgetl(fid)); end
first=sscanf(zeile,'%f');
mm=fscanf(fid,'%f',[length(first),Inf])';
fclose(fid);
mm(end+1,1:length(first))=first(:)';
if length(first)<5, % xm zm rho
    zm=unique(mm(:,2));
    z=0;
    for i=1:length(zm),
        z=[z z(i)+2*(zm(i)-z(i))];
    end
    xm=unique(mm(:,1));
    x=1.5*xm(1)-0.5*xm(2);
    for i=1:length(xm),
        x=[x x(i)+2*(xm(i)-x(i))];
    end
    M=zeros(length(xm),length(zm));
    [xx,ix]=ismember(mm(:,1),xm);
    [zz,iz]=ismember(mm(:,2),zm);
    nr=3;
else
    x=unique(mm(:,1));x(end+1)=max(mm(:,2));    
    z=unique(mm(:,3));z(end+1)=max(mm(:,4));    
    M=zeros(length(x)-1,length(z)-1);
    [xx,ix]=ismember(mm(:,1),x);
    [zz,iz]=ismember(mm(:,3),z);
    nr=5;
end
for l=1:min(length(ix),length(iz)),
    M(ix(l),iz(l))=mm(l,nr);    
end
Cov=[];
if ((nargout>3)|(nargout==1))&&(size(mm,2)>nr), %coverage
    Cov=M;
    for l=1:min(length(ix),length(iz)),
        Cov(ix(l),iz(l))=mm(l,nr+1);    
    end
end
if nargout==1, % model structure
    MM=M;Lay=0;
    M=struct('M',MM,'x',x,'z',z,'Lay',Lay,'Cov',Cov,'R',[]);    
end   

function zeile=destrip(zeile)
% strip string from comments (with # character)
aa=strfind(zeile,'#');
if ~isempty(aa), zeile=zeile(1:aa(1)-1); end
