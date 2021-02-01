function showsipdata1d(ab2,f,A,ind,cax)

% SHOWSIPDATA1D - Show SIP 1d data as image plot
% showsipdata1d(ab2,f,A,toshow,caxis)

if (nargin<4)||isempty(ind), ind=logical(ones(size(A))); end
if nargin<5, cax=minmax(A(ind(:))); end
if prod(cax)<0, % diff. signs
    cax=[-1 1]*max(abs(cax));
%     cmap=b2r(64);
    cmap=jet(64);
else
    cmap=jet(64);
end
cind=getcindex(A(:),cax(1),cax(2),size(cmap,1));
cla;caxis(cax);
nr=size(A,1);
nf=size(A,2);
for i=1:nr,
    for j=1:nf,
        if ind(i,j),
            patch([-1 -1 1 1]*0.5+j,[-1 1 1 -1]*0.5+i,...
                cmap(cind(nr*(j-1)+i),:));
        end
    end
end
xlim([0.5 nf+0.5]);
ylim([0.5 nr+0.5]);
% imagesc(deltaPhi);
% alpha(double(P>0));
cb=colorbar;
if cax(1)==-cax(2), 
    cba=get(cb,'Children');
%     set(cba,'CData',cmap); 
end
axis ij;
set(gca,'YTick',1:length(ab2),'YTickLabel',num2strcell(ab2),...
    'XTick',1:length(f),'XTickLabel',strrep(num2strcell(rndig(f)),'000','k'));
xlabel('f in Hz');
ylabel('AB/2 in m');
title('\phi in mrad');
