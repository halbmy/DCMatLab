[mids,seps,ii,kk,midpoint,dm]=midkonf2d(N);
xx=repmat(midpoint,1,4)+repmat([-1 1 1 -1]*dm/2,length(midpoint),1);
yy=repmat(kk,1,4)+repmat([-1 -1 1 1]/2,length(kk),1);
clf;cmap=jet(64);
aa=log10(N.r);mi=min(aa);ma=max(aa);
jj=round((aa-mi)/(ma-mi)*(size(cmap,1)-1))+1;
jj(jj<1)=1;jj(jj>size(cmap,1))=size(cmap,1);
% p=patch(xx,yy,aa,'Faces',cmap(jj,:))
clf;for i=1:length(jj),
    patch(xx(i,:),yy(i,:),cmap(jj(i),:),'Linestyle','none','UserData',i);
end
set(gca,'YDir','reverse');axis tight

vor='pppddpweslddddgr';
yt=get(gca,'YTick');
if find(yt-round(yt)), %0.5 etc
    yt=round(min(yt)):floor(max(yt));
    set(gca,'YTick',yt);
end
for i=1:length(yt),
    yy=seps(yt(i));
    switch yy,
        case 9999, % circulated dipole
            ytli='cc';
        case 30001, % wenner=schlumberger sep.1
            ytli='ws1';
        otherwise
            st=fix(yy/10000)*2+1;
            if st<length(vor)-1, ytli=vor(st:st+1); else ytli='gr'; end
            aa=mod(yy,10000);
            bb=fix(aa/100+1);
            if bb>1, ytli=[ytli num2str(bb) '-']; end
            if st==15, ytli(end)='+'; end
            %                 if st==17, cc=mod(100-aa,100); else 
            cc=mod(aa,100); 
            %                 end
            ytli=[ytli num2str(cc)];
        end
        ytl(i,1:length(ytli))=ytli;
    end
    set(gca,'YTickMode','manual','YTickLabel',ytl);