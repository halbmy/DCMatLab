function showpseudoflat(N,feld,cmin,cmax,iscbar,islog)

% SHOWPSEUDOFLAT - Show flat pseudosection
% showpseudoflat(N,field,cmin,cmax,iscbar,islog)

cla reset;
if nargin<2, feld=N.r; end
if nargin<5, iscbar=1; end
if nargin<6, islog=(min(feld)>0); end
if ischar(N), N=readinv2dfile(N); end
if nargin<2, feld=N.r; end
if nargin<3, cmin=min(feld); end
if nargin<4, cmax=max(feld); end
if islog, feld=log10(feld);cmin=log10(cmin);cmax=log10(cmax); end
nel=size(N.elec,1);
aa=N.a-1;bb=N.b-1;mm=N.m-1;nn=N.n-1;
fi=find(bb<aa);bb(fi)=bb(fi)+nel;
fi=find(mm<bb);mm(fi)=mm(fi)+nel;
fi=find(nn<mm);nn(fi)=nn(fi)+nel;
sep=abs(mm-bb);
sep=min(sep,nel-2-sep);
mid=mod(bb+sep/2,nel)+1;
fi=find(mm-bb>(nel-2)/2);
mid(fi)=mod(aa(fi)-sep(fi)/2,nel)+1;
cmap=colormap(jet);lcm=length(cmap)-1;
xx=[1 1 -1 -1 1]/2;yy=[1 -1 -1 1 1]/2;
for i=1:length(feld),
    if isfinite(feld(i)),
        cind=round(1+(feld(i)-cmin)/(cmax-cmin)*lcm);
        if cind<1, cind=1; end
        if cind>lcm, cind=lcm; end
        pa=patch(mid(i)+xx,sep(i)+yy,cmap(cind,:),'EdgeColor','black');%cmap(cind,:));
        if mid(i)<5, pa=patch(mid(i)+nel+xx,sep(i)+yy,cmap(cind,:),'EdgeColor','black'); end
    end
end
line([1 1]*(nel+0.5),[min(sep)-0.5 max(sep)+0.5],'Color','red');
axis equal tight
set(gca,'YDir','reverse','XAxisLocation','top');
xlabel('Electrode position');
ylabel('Separation factor');
caxis([cmin cmax]);
if iscbar,
    hc=colorbar('horiz');
    xt=get(hc,'XTick');xtl=num2strcell(xt);
    if islog,
        if max(feld(:))>2, % >100
            xtl=num2strcell(round(10.^xt));
        elseif max(feld(:))>1, % >10
            xtl=num2strcell(round(10.^xt*10)/10);
        else
            xtl=num2strcell(round(10.^xt*100)/100);
        end
        set(hc,'XTickMode','manual','XTickLabel',xtl);
    end
end


function ce=num2strcell(vec)

for i=1:length(vec),
    ce{i}=num2str(vec(i));
end