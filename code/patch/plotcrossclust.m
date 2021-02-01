function plotcrossclust(data,cl,vv)

cols={'br','bgr','bcgy','bcgyr','bcgyrm','bcgyrmk'};
clf;hold on;
maxcl=max(cl);
for j=1:maxcl,
    fi=find(cl==j);
    if ~isempty(fi),
        plot(data(fi,1),data(fi,2),[cols{maxcl-1}(j) '.']);
    end
    hold on;
end
if nargin>2,    
    xl=xlim;yl=ylim;
    for j=1:size(vv,1), 
        plot(vv(j,1),vv(j,2),[cols{maxcl-1}(j) 'x'],'MarkerSize',10);
        plot(vv(j,1),yl(1),[cols{maxcl-1}(j) '+'],'MarkerSize',10);
        plot(xl(1),vv(j,2),[cols{maxcl-1}(j) '+'],'MarkerSize',10);
    end
end
% for i=1:size(data,1),
%     plot(data(:,1),data(:,2),[cols{end}(cl(i)) '.']);
% end
axis tight;grid on;hold off;