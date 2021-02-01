function Bor = loadborfile(borfile,z)

% LOADBORFILE - Load borehole (stratigraphic) file
% Bor = loadborfile(filename)
% Bor..struct containing
%    pos - borehole position vector
%    lay - cell array of layer vectors

fid=fopen(borfile,'r');
if fid<0, error('File does not exist!'); end
n=str2double(destrip(fgetl(fid)));
Bor.pos=zeros(n,1);
for i=1:n,
    aa=str2num(destrip(fgetl(fid)));
    Bor.lay{i}=[];
    Bor.pos(i)=aa(1);
    for j=1:aa(2),
        bb=str2double(destrip(fgetl(fid)));
        if bb>0, Bor.lay{i}(end+1)=bb; end
    end
end
fclose(fid);
if nargin>1, % z given -> calculate nlay
    for i=1:length(Bor.lay),
        Bor.nlay{i}=zeros(size(Bor.lay{i}));
        for j=1:length(Bor.lay{i}),
            [mi,k]=min(abs(z-Bor.lay{i}(j)));
            Bor.nlay{i}(j)=k(1);
        end
        Bor.nlay{i}(diff(Bor.nlay{i})==0)=[];
    end
end
return
dx=10;
for i=1:length(Bor.lay),
    for j=1:length(Bor.lay{i}),
        line([-1 1]*dx+Bor.pos(i),Bor.lay{i}(j)*[1 1],'Color','w');
        %       line([-1 1]*dx*0.8+Bor.pos(i),z(Bor.nlay{i}(j))*[1 1],'Color','b');
    end
end
return
% function zeile=destrip(zeile)
% % strip string from comments (with # character)
% aa=strfind(zeile,'#');
% if ~isempty(aa), zeile=zeile(1:aa(1)-1); end
