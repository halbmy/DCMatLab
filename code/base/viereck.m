function drawrectangle(x1,x2,z1,z2,fillc)

% DRAWRECTANGLE - draw (filled) rectangle
% drawrectangle(x1,x2,z1,z2[,fillcolor])

if nargin<5, fillc=''; end

lix=[x1 x2 x2 x1 x1];
liy=[z1 z1 z2 z2 z1];
if isempty(fillc),
    line(lix,liy,'LineWidth',1,'Color','black');
else
    patch(lix,liy,fillc,'LineWidth',3);
end