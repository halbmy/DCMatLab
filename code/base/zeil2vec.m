function [id,vec]=zeil2vec(zeile)

% ZEIL2VEC - reads line from string as follows
% identifier=number1/number2/...
% [identifier,vector_of_numbers]=zeil2vec(string)

id=zeile;vec=0;
ig=strfind(zeile,'=');
if isempty(ig), return; end
id=zeile(1:ig-1);
zeile=zeile(ig+1:end);
vec=[];
while strfind(zeile,'/'),
    pos=strfind(zeile,'/');
    aa=str2num(strrep(zeile(1:pos(1)-1),',','.'));
    if isempty(aa), vec(end+1)=0; else vec(end+1)=aa; end
    if pos>=length(zeile), return; end
    zeile=zeile(pos(1)+1:end);
end
if ~isempty(zeile), 
    aa=str2num(strrep(zeile,',','.')); 
    if ~isempty(aa), vec(end+1)=aa; end
end