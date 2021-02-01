function plotholes(holes,dx2,cols)

% PLOTHOLES - Plot borehole information over plot
% plotholes(holeinfo,width,colors)
% holeinfo - matrix of rows for each hole containing
%            position height depth1 depth2 ...
% width - width of boreholes to be plotted
% colors - cell of one color characters {'y','b','r'}

if nargin==0, return; end
if nargin<2, dx2=2; end
if nargin<3, cols={'y','b','r','c','g'}; end

if ischar(cols),%
   old=cols;
   cols={};for i=1:length(old), cols{i}=old(i); end
end
% yl=ylim;
for i=1:size(holes,1),
    xx=holes(i,1)+[-1 1 1 -1]*dx2;
    set(patch(xx,holes(i,2)-[0 0 1 1]*holes(i,3),cols{1}),'EraseMode','None');
%     set(patch(xx,holes(i,2)-holes(i,3)-[0 0 1 1]*(holes(i,4)-holes(i,3)),cols{2}),'EraseMode','None');
    %         set(xx,[(holes(i,2)-holes(i,4))*[1 1] yl(1)*[1 1]],cols{3}),'EraseMode','None');
    for j=3:7,
        if size(holes,2)>j,
            set(patch(xx,[(holes(i,2)-holes(i,j))*[1 1] (holes(i,2)-holes(i,j+1))*[1 1]],cols{j-1}),'EraseMode','None');
        end
    end

%     if size(holes,2)>4,
%         set(patch(xx,[(holes(i,2)-holes(i,4))*[1 1] (holes(i,2)-holes(i,5))*[1 1]],cols{3}),'EraseMode','None');
%     end
%     if size(holes,2)>5,
%         set(patch(xx,[(holes(i,2)-holes(i,5))*[1 1] (holes(i,2)-holes(i,6))*[1 1]],cols{4}),'EraseMode','None');
%     end
%     if size(holes,2)>5,
%         set(patch(xx,[(holes(i,2)-holes(i,5))*[1 1] (holes(i,2)-holes(i,6))*[1 1]],cols{4}),'EraseMode','None');
%     end
end