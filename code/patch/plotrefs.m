function plotrefs(filename,color)

% PLOTREFS - Plot refractor lines
% plotrefs(filename)

% filename='lay.xy';

if nargin<1, error('Specify filename!'); end
if nargin<2, color='white'; end
fid=fopen(filename);
A=[];
hold on;
while 1,
    zeile=fgetl(fid);
    if ~ischar(zeile), 
        if ~isempty(A), line(A(:,1),A(:,2),'Color',color); end
        break; 
    end
    po=strfind(zeile,'#');
    if ~isempty(po), zeile=zeile(1:po(1)-1); end
    if isempty(zeile)
        if ~isempty(A), line(A(:,1),A(:,2),'Color',color); end
        A=[];
    else
        xz=str2num(zeile);
        if length(xz)>1, 
            A(end+1,1:2)=xz(1:2); 
        else
            line(A(:,1),A(:,2),'Color',color);
            A=[];
        end
    end
end
hold  off;
fclose(fid);
