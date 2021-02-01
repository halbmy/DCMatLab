function Shot = readpof(filename)

% READPOF - Read promax output file (POF)
% Shot = readpof(filename)

if nargin<1, error('Specify filename'); end
%filename='c:\halbmy\2d\sollstedt\fb_2.pof';
% A=textread(filename,'','headerlines',5);
Shot=[];
l=0;sh=[];ge={};Shot.t=[];
fid=fopen(filename,'r');
if fid<0, error(['File ' filename ' could not be opened']); end
zeile=fgetl(fid);
while ischar(zeile),
    aa=str2num(zeile);
    if length(aa)==3, % new shot
        l=l+1;
        sh(l)=aa(1);
        ge{l}=aa(2);
        Shot.loc(l)=aa(1);
        Shot.x{l}=aa(2);
        Shot.tt{l}=aa(3);
        Shot.t(end+1)=aa(3)/1000;
    end
    if (l>0)&&(length(aa)==2), % new position
        ge{l}(end+1)=aa(1);
        Shot.x{l}(end+1)=aa(1);
        Shot.tt{l}(end+1)=aa(2);
        Shot.t(end+1)=aa(2)/1000;
    end
    zeile=fgetl(fid);
end
fclose(fid);
pos=sh(:);Shot.t=Shot.t(:);
for i=1:length(ge), pos=[pos;ge{i}(:)]; end
Shot.pos=unique(pos);
for i=1:length(ge), 
    Shot.ns{i}=find(Shot.pos==sh(i));
    [tf,Shot.nx{i}]=ismember(ge{i},Shot.pos);
end
Shot.pos(:,2)=0;Shot.locz=zeros(size(Shot.loc));
writetom(Shot,'fb2.tom');