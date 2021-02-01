function Shot=readunishot(filename)

% READUNISHOT - Read shot from unified data format
% Shot = readunishot(filename)

N=readunifile(filename);
if isfield(N,'x'), 
    Shot.pos=N.x;
    if isfield(N,'y'), Shot.pos(:,2)=N.y; end
    if isfield(N,'z'),
        Shot.pos(:,end+1)=N.z;
    else
        if isfield(N,'h'), Shot.z=N.h; end
        if isfield(Shot,'z')&&isfield(N,'d')&&(length(Shot.z)==length(N.d)),
            Shot.z=Shot.z-N.d; end
    end
else
    Shot.pos=N.elec;
end
if ~isfield(N,'s')&isfield(N,'a'), N.s=N.a; end
if ~isfield(N,'g')&isfield(N,'m'), N.g=N.m; end
us=unique(N.s);
tf=isfield(N,'t');
if tf, Shot.t=N.t; end
if isfield(N,'g'), Shot.g=N.g; end
if isfield(N,'s'), Shot.s=N.s; end
for i=1:length(us),
    Shot.ns{i}=us(i);
    fi=find(N.s==us(i));    
%     [Shot.nx{i},I,J]=unique(N.g(fi));
    [Shot.nx{i},I]=sort(N.g(fi));
    if tf, Shot.tt{i}=Shot.t(fi(I))*1000; end
    Shot.nn{i}=fi(I);
end
% if tf, 
%     Shot.t=[];
%     for i=1:length(Shot.tt), Shot.t=[Shot.t;Shot.tt{i}]; end
%     Shot.t=Shot.t/1000;
% end
