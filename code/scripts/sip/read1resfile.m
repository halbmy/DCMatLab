function [f,rhoa,phi,drhoa,dphi] = read1resfile(fname,readsecond)

% READ1RESFILE - Read single RES (Radic SIP) file
% Data = read1resfile(fname)
% [f,rhoa,phi,drhoa,dphi] = read1resfile(fname)

fid=fopen(fname);
f=[];phi=[];dphi=[];rhoa=[];drhoa=[];
zeile=strrep(fgetl(fid),'"','');
while ischar(zeile),
    if (length(zeile)>3)&&isequal(zeile(1:4),'Freq'), break; end
    zeile=strrep(fgetl(fid),'"','');
end
if (nargin>1)&&(readsecond),
    zeile=strrep(fgetl(fid),'"','');
    while ischar(zeile),
        if (length(zeile)>3)&&isequal(zeile(1:4),'Freq'), break; end
        zeile=strrep(fgetl(fid),'"','');
    end
end
zeile=fgetl(fid);
while ischar(zeile)&&(~isempty(zeile)),
    zeile=strrep(zeile,'"','');
    aa=str2num(zeile);
    f(end+1)=aa(1);
    rhoa(end+1)=aa(2);%*k;
    phi(end+1)=aa(3);
    drhoa(end+1)=aa(4);
    dphi(end+1)=aa(5);
    zeile=fgetl(fid);
end
fclose(fid);
phi=-phi*pi/180;dphi=dphi*pi/180;
[f,so]=sort(f);
rhoa=rhoa(so);
phi=phi(so);
drhoa=drhoa(so);
dphi=dphi(so);
if nargout==1,
   ff=f;
   f=struct('f',ff,'rhoa',rhoa,'phi',phi,'drhoa',drhoa,'dphi',dphi);
end
