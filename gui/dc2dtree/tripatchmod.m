function [cmin,cmax]=tripatchmod(Mesh,att,N,cmin,cmax,canot)

%  tripatchmod(Mesh,att,N);

if nargin==0, error('No mesh specified!'); end
if nargin<2, % no att specified
    patch('Vertices',Mesh.node,'Faces',Mesh.cell,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]);
    axis equal tight
    return
end
if (nargin<2)||(length(att)~=Mesh.ncells), att=Mesh.cellattr; end
%     rand(Mesh.ncells,1)*2-1; end
if nargin<4, cmin=1; end
if (nargin==4)&&ischar(cmin), canot=cmin;cmin=1; end
if (nargin==3)&&ischar(cmin), canot=N; end
if nargin<5, cmax=1; end
if nargin<6, canot=0; end
nel=size(Mesh.cell,1);
if length(att)<nel, att(nel)=1; end
if length(att)>nel, att(nel+1:end)=[]; end
islog=(min(att)>0)&&(cmin>0);
if islog, 
    att=log10(att); 
    cmin=log10(cmin);
    cmax=log10(cmax);
end
if nargin<4, 
    [NN,X]=hist(att,100);
    C=cumsum(NN)/sum(NN);
    cmin=X(min(find(C>0.05)));%min(att);
    if isempty(cmin), cmin=min(att); end
end
if nargin<5,
    [NN,X]=hist(att,100);
    C=cumsum(NN)/sum(NN);
    cmax=X(max(find(C<0.95)));%max(att);
    if isempty(cmax), cmax=max(att); end
end
if ~(cmax>cmin), cmax=cmin+1; end
if islog|(max(att)*min(att)>0),
    cmap=colormap(jet);
else
    cmap=colormap(b2r);
    cmax=max(abs([cmin cmax]));cmin=-cmax;
end
lcm=length(cmap)-1;
emap=cmap;
if (length(att)<100)|(length(unique(att))<10), emap(:)=0; end
cla reset;
if 1,
for i=1:length(att),
    cind=round(1+(att(i)-cmin)/(cmax-cmin)*lcm);
    if cind<1, cind=1; end
    if cind>lcm, cind=lcm; end
%     patch(NODE(ELE(i,2:4),1),NODE(ELE(i,2:4),2),cmap(cind,:),'EdgeColor',emap(cind,:));
    patch(Mesh.node(Mesh.cell(i,1:3),1),Mesh.node(Mesh.cell(i,1:3),2),cmap(cind,:),'EdgeColor',emap(cind,:));
end
else
    edgecolor='none';if length(unique(att))<10, edgecolor='black'; end
    patch('Vertices',Mesh.node,'Faces',Mesh.cell,'FaceVertexCData',...
       att(:),'FaceColor','flat','EdgeColor',edgecolor);%faces,'EdgeColor', edges );    
end
axis equal tight
caxis([cmin cmax]);
cb=colorbar('horiz');
xtl=cellstr(get(gca,'XTickLabel'));
xtl{end-1}='x/m';
set(gca,'XTickLabel',xtl);
ytl=cellstr(get(gca,'YTickLabel'));
ytl{end-1}='y/m';
set(gca,'YTickLabelMode','manual','YTickLabel',ytl);
% set(cb,'DataAspectRatio',[1 32*2 1]);
if islog,
    xt=get(cb,'XTick');
    if (cmax-xt(end))/(cmax-cmin)>0.05,  %enough space left
        xt(end+1)=cmax; 
        set(cb,'XTick',xt);
    end
    xt=rndig(10.^xt);
    xtl=num2strcell(xt);
    if ischar(canot), xtl{end-1}='Ohmm'; end
    set(cb,'XTickLabelMode','manual','XTickLabel',xtl);
end
dar=get(cb,'DataAspectRatio');
set(cb,'DataAspectRatio',dar.*[1 32 1]);
if (nargin>2)&&isfield(N,'elec'),
    elec=[N.elec];%;N.elec(1,:)];
    hold on;plot(elec(:,1),elec(:,2),'kx-');
    plot(elec(1,1),elec(1,2),'ro');hold off
else
    rad=max(abs([xlim ylim]));
    if rad<=1, %obviously a tree
        nel=24;
        elec=rad*sin((0:nel)'/nel*2*pi);
        elec(:,2)=rad*cos((0:nel)'/nel*2*pi);
    end
%     hold on;plot(elec(:,1),elec(:,2),'kx-');
%     plot(elec(1,1),elec(1,2),'ro');hold off
end
if islog, cmin=10^cmin;cmax=10^cmax; end