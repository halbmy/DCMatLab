function iconify(h,ico)

% ICONIFY - Put Icon on figure
% iconify(h,iconfile)
% h - figure handle
% iconfile - icon file

if ~ispc, return; end
if nargin<2, ico='.\dc2dtree.ico'; end
if ishandle(h)&&exist(ico,'file'),
    name=get(h,'Name');
    if ~isempty(name), icon(101,name,ico); end
end