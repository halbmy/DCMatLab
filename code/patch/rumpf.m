fprintf('\nLoaded %d values, min=%g max=%g\n',length(att),min(att),max(att));
tetdef=[1 2 3;1 2 4;1 3 4;2 3 4];
[x,y,z]=en2mid(ELE,NODE);
if find([sx sy sz]), % Schnitte
    clog=tetvis(x,y,z,att,[1 1 0.5],[],sx,sy,sz);
else % testing reasons
    clf;plot3(0,0,0);set(gca,'ZDir','normal');%reverse
    clog=(min(att)>0);
    caxis(log10([100 500]));
end
% xx=NODE(:,1);yy=NODE(:,2);zz=NODE(:,3);
xx=x;yy=y;zz=z;
xlim([min(xx) max(xx)]);
ylim([min(yy) max(yy)]);
zlim([min(zz) max(zz)]);
cmap=colormap;cach=caxis;
%if clog, cach=10.^cach; end
if drawfaces, %draw faces
    ca=exist([meshname '.e']); % carsten-format (ohne 1. Zeile und Spalte)
    if ca,
        %fid=fopen([meshname '.f'],'r');
        %FACE=fscanf(fid,'%f',[6 Inf])';
        %fclose(fid);
	FACE=load([meshname '.f'],'-ascii');
        FACE(:,4:5)=[];
    else
        fid=fopen([meshname '.face'],'r');
        si=fscanf(fid,'%d',[2 1]);
        FACE=fscanf(fid,'%f',[5 si(1)])';
        fclose(fid);
        FACE(:,1)=[];
    end
    fi=find(FACE(:,4)==-1); 
    %==-1 for neumann, ==-2 for dirichlet, <0 for both
    FACE=FACE(fi,1:3)+1;
    rf=[1 2 3 1];
    hold on
    for e=1:size(FACE,1),
        pts=NODE(FACE(e,:),1:3);
        plot3(pts(rf,1),pts(rf,2),pts(rf,3),'k');
    end
    hold off
end
for l=1:length(isoval),
    gtlt='less';
    if isoval(l)>0,
        ind=find(att>isoval(l));gtlt='greater';
    else
        ind=find(att<-isoval(l));
    end
    fprintf('Found %d values %s than %g\n',length(ind),gtlt,abs(isoval(l)));
    for i=1:length(ind),
        pts=NODE(ELE(ind(i),1:4),1:3);
        col=[1 0 0];val=att(ind(i));if clog, val=log10(val); end
        cind=1+round((val-cach(1))/(cach(2)-cach(1))*(length(cmap)-1));
        if cind>length(cmap), cind=length(cmap); end
        if cind<1, cind=1; end
        col=cmap(cind,:);
        for j=1:4,
            td=tetdef(j,:);
            patch(pts(td,1),pts(td,2),pts(td,3),col);
        end
    end
end
xy='Y';cbdir='vert';
cbdir='horiz';
if isequal(cbdir,'horiz'), xy='X'; end 
cb=colorbar(cbdir);
if clog,
   yt=10.^(get(cb,[xy 'Tick']));
   yt=round(yt);
   ytl=num2str(yt(:));
%    ytl(end-1,1:
    set(cb,[xy 'TickLabel'],ytl);
end
daspect([1 1 1/2]);
zlim([0 6]);
xlabel('x in m');
ylabel('y in m');
zlabel('z in m');
view(3);
