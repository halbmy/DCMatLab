function Shot = readmax(filename)

% READMAX - Read refraction file (Max format)
% Shot = readmax(filename)

if nargin<1, error('Specify filename!'); end
% filename='max.dat';
fid=fopen(filename,'r');
head1=fgetl(fid);
head2=fgetl(fid);
fclose(fid);
i=0;cols=[];pos=[];
rest=head2;
while ~isempty(rest),
    i=i+1;
    [tok,rest]=strtok(rest);
    if any(tok)&&isequal(lower(tok(1:2)),'sp'),
        cols(end+1)=i;
        pos(end+1)=str2num(tok(3:end));
    end
end
A=textread(filename,'','headerlines',2);
gpos=A(:,2);
if size(A,2)>max(cols), % topo present
    gpos(:,2)=A(:,max(cols)+1);fi=find(gpos(:,2));fi1=find(gpos(:,2)==0);
    if any(fi)&any(fi1), gpos(fi1,2)=interp1(gpos(fi,1),gpos(fi,2),gpos(fi1,1),'linear','extrap'); end
else gpos(:,2)=0; end
posz=interp1(gpos(:,1),gpos(:,2),pos,'linear','extrap');
Shot.pos=unique([gpos;pos(:) posz(:)],'rows');
[tf,loc]=ismember(gpos(:,1),Shot.pos(:,1));
Shot.t=[];Shot.loc=pos(:);Shot.locz=posz(:);
% [tf,Shot.ns]=ismember(pos,Shot.pos(:,1));
for i=1:length(cols),
    [tf,Shot.ns{i}]=ismember(pos(i),Shot.pos(:,1));
    fi=find(isfinite(A(:,cols(i))));
    Shot.nx{i}=loc(fi);
    Shot.t=[Shot.t;A(fi,cols(i))/1000];
    Shot.x{i}=gpos(fi,1);
    Shot.tt{i}=A(fi,cols(i));
end
Shot.pos=round(Shot.pos*100)/100;