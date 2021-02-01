function plotshot(Shot,tvgl,offset)

% PLOTSHOT - Plot seismic shots and traveltimes
% plotshot(Shot); %with Shot is structure
% plotshot(Shot,t_compare)

if nargin<3, offset=[]; end
cla;
cols='bgrcmyk';
set(gca,'YTickMode','auto','YTickLabelMode','auto');
l=0;m1='x-';
if length(Shot.t)>5000, m1='.-'; end
if nargin>1, 
    m1='x';m2='-'; 
    if length(Shot.t)>5000, m1='.'; end
end
% if ~isfield(Shot.loc), Shot.loc
if isfield(Shot,'s')&&isfield(Shot,'g'), % only s and g
    unis=unique(Shot.s);
    if ~isfield(Shot,'pos')&&isfield(Shot,'elec'), Shot.pos=Shot.elec; end
    for i=1:length(unis),
       col=cols(mod(i-1,7)+1);
       fi=find(Shot.s==unis(i));
       plot(Shot.pos(Shot.g(fi),1),Shot.t(fi)*1000,[col m1]);
       if nargin>1,
           plot(Shot.pos(Shot.g(fi),1),tvgl(fi)*1000,[col m2]);
       end
       if i==1, hold on; end
       to=0;
       if length(offset)>=i, to=offset(i)*1000; end
       plot(Shot.pos(unis(i),1),to,[col '*']);
    end
elseif isfield(Shot,'ns')&&isfield(Shot,'nx'),
    for i=1:length(Shot.ns),
        col=cols(mod(i-1,7)+1);
        %    plot(Shot.x{i},Shot.tt{i},[col m1]);
        plot(Shot.pos(Shot.nx{i},1),Shot.tt{i},[col m1]);
        if i==1, hold on; end
        %    plot(Shot.loc(i),0,[col '*']);
        to=0;
        if length(offset)>=i, to=offset(i)*1000; end
        plot(Shot.pos(Shot.ns{i},1),to,[col '*']);
        %just to check shot numbers
        %    text(Shot.pos(Shot.ns{i},1),to,num2str(i),'FontSize',8,'Color',col);
        if nargin>1,
%             lx=length(Shot.nx{i});            
            xx=Shot.pos(Shot.nx{i},1);
            tt=tvgl(l+1:l+lx);
            plot(xx,tt,[col m2]);
%             l=l+lx;
        end
    end
end
hold off;
set(gca,'xlim',[min(Shot.pos(:,1)) max(Shot.pos(:,1))]);
% xlabel('shot/receiver position [m]');ylabel('travel time [ms]');
grid on;box on;
set(gca,'YDir','reverse','XAxisLocation','top');
xtl=cellstr(get(gca,'XTickLabel'));xtl{end-1}='x/m';
ytl=cellstr(get(gca,'YTickLabel'));ytl{end-1}='t/ms';
set(gca,'XTickLabel',xtl,'YTickLabel',ytl);