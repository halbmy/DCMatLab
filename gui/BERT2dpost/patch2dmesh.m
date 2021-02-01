function patch2dmesh(Mesh,field,alphavalue)

% PATCH2DMESH - Patch a 2D mesh
% patch2dmesh(Mesh)
% patch2dmesh(Mesh,field)
% patch2dmesh(Mesh,field,alphavalues)

if ~isfield(Mesh,'dim')||(Mesh.dim~=2)||~isfield(Mesh,'node')||~isfield(Mesh,'cell'),
    error('Mesh is not a valid 2d mesh!'); end
if ~isfield(Mesh,'nnodes'), Mesh.nnodes=size(Mesh.node,1); end
if ~isfield(Mesh,'ncells'), Mesh.ncells=size(Mesh.cell,1); end
if nargin<2,
    if isfield(Mesh,'ZZcellattr'), field=Mesh.cellattr; else
        field=(1:Mesh.ncells)'; end
end
if nargin<3, alphavalue=ones(size(field)); end
if length(field)>Mesh.ncells, error('Field too long!'); end
islog=0;
if min(field)>0, field=log10(field);islog=1; end
cmin=min(field);cmax=max(field);
perc=5;
[N,X]=hist(field,100);
C=cumsum(N)/sum(N);
% cmin=X(min(find(C>perc/100)));
% cmax=X(max(find(C<1-perc/100)));

cmap=colormap;emap=cmap;lecm=length(cmap)-1;
for i=1:length(field),
    cind=round(1+(field(i)-cmin)/(cmax-cmin)*lecm);
    if cind<1, cind=1; end
    if cind>lecm, cind=lecm; end
    col=cmap(cind,:);
    col=col*alphavalue(i)+1-alphavalue(i);
    nums=Mesh.cell(i,:);
    patch(Mesh.node(nums,1),Mesh.node(nums,2),col,'EdgeColor',...
	emap(cind,:),'LineStyle','none');
        %[0 0 0]);
end
axis equal tight
xtl=cellstr(get(gca,'XTickLabel'));xtl{end-1}='x/m';
ytl=cellstr(get(gca,'YTickLabel'));ytl{end-1}='z/m';
set(gca,'XTickLabel',xtl,'YTickLabel',ytl);

cb=colorbar('horiz');
% xt=get(cb,'XTick');
xt=linspace(cmin,cmax,10);
if islog, xt=10.^xt; end
fi=find(abs(xt)>1);xt(fi)=round(xt(fi)*10)/10;
fi=find(abs(xt)>10);xt(fi)=round(xt(fi));
set(cb,'XTick',linspace(0,1,10),'XTickLabel',num2strcell(xt));
