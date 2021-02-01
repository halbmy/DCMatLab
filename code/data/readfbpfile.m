function Shot = readfbpfile(filename,dx)
% READFBPFILE - Read FBP (first break picks) file
% Shot = readfbpfile(filename)

% filename='wuelf.fbp';
if nargin<2, dx=0; end
pp=fileparts(filename);
fid=fopen(filename,'r');
maxg=0;i=0;pos=[];
zeile=fgetl(fid);
while isstr(zeile),
    i=i+1;
    %%
    trfile=sscanf(zeile,'%s%*s%*s');
    if isempty(strfind(trfile,'.pick')), trfile=[trfile '.pick']; end
    pos(i)=sscanf(zeile,'%*s%f%*s');
    add=sscanf(zeile,'%*s%*s%d');
    if isempty(add), add=0; end
    %%
    A=load(fullfile(pp,trfile));
    A(:,2)=A(:,2)+add;
    maxg=max(maxg,max(A(:,2)));
    A=flipud(A);
    [U,I,J]=unique(A(:,2));%,'last');
    TN{i}=A(I,:);
    zeile=fgetl(fid);
end
%%
fclose(fid);
%%
if dx==0, dx=round((max(pos)-min(pos))/maxg); end
allpos=[(0:maxg-1)'*dx;pos(:)];
Shot=[];Shot.t=[];
[Shot.pos,I,J]=unique(allpos);
Shot.pos(:,2)=0;
for i=1:length(pos),
    Shot.ns{i}=J(maxg+i);
    Shot.nx{i}=J(TN{i}(:,2));
    Shot.tt{i}=TN{i}(:,1)*1000;
    Shot.t=[Shot.t;Shot.tt{i}/1000];
end