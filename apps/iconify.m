function h=iconify(h,ico)

% ICONIFY - Put Icon on figure
% iconify(h,iconfile)
% h - figure handle
% iconfile - icon file

return
if nargin<2, ico='dc2dinvres.ico'; end
if ispc&&ishandle(h)&&(exist(ico)==2),
    name=get(h,'Name');
    if ~isempty(name), icon(101,name,ico); end
%     if ~isempty(name), icon(101,name,which(ico)); end
end