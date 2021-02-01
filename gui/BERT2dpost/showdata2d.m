function [mids,seps,ii,kk]=showdata2d(NN,feld,MAL)

% SHOWDATA2D show datum points
% showdata2d(N,[,field[,DRAWOPTS])
% N..Structure of electrode numbers(a,b,m,n), 
% k-factors(k), measurements(r) and
% elec .. electrode positions
% field..data to plot, otherwise N.r

if nargin<1, error('Data structure as input argument required!'); end
if nargin<2, feld=NN.r; end
if isempty(feld), feld=NN.r; end
if nargin<3,
  MAL=struct('cauto',1,'cmap',0);      
end
clog=1;
if ~isfield(MAL,'cauto')&isfield(MAL,'cmin')&isfield(MAL,'cmax'), MAL.cauto=0; end
if isfield(MAL,'clog'), MAL.log=MAL.clog; end
if isfield(MAL,'xdir'), xdir=MAL.xdir; else xdir=0; end
if isfield(MAL,'cmax'), cmax=MAL.cmax; else cmax=1; end
if isfield(MAL,'cmin'), cmin=MAL.cmin; else cmin=0; end
if isfield(MAL,'log'), clog=MAL.log; end
if min(feld(:))<0, clog=0; end
if isfield(MAL,'cmap'), 
    cmap=MAL.cmap; 
else 
    cmap=((min(feld)<0)&(clog==0))*2; 
end
if isfield(MAL,'cauto'), cauto=MAL.cauto; else cauto=1; end
if isfield(MAL,'cbar'), cbar=MAL.cbar; else cbar=1; end
if length(unique(NN.elec(:,2)))==1, NN.elec(:,2)=0; end
if 0&&length(find(NN.elec(:,2)>0))>size(NN.elec,1)/2, %topo in elec
    % check for tape correction
    NN.elec(:,2)=0;
end
if find(NN.elec(:,2)~=0), % borehole electrodes
    sel=find(NN.elec(:,2)==0);
    if isempty(sel), % ist mir unklar
%         N=[];N.a=[];
        N=NN;
        N.elec(:,2)=0; 
    else
        fi=findmess(NN,sel,sel);
        N.elec=NN.elec;
        N.a=NN.a(fi);N.b=NN.b(fi);
        N.m=NN.m(fi);N.n=NN.n(fi);
        N.k=NN.k(fi);N.r=NN.r(fi);
        field=feld;
        feld=field(fi);
    end
else % surface measurement
    N=NN;
end
mids=[];seps=[];
if length(N.a)>0,
%     if isfield(N,'eind'),
%         for l=1:length(N.eind),
%             si=size(N.eind{l},1);
%             data(1:si,l)=log10(N.eind{l}(:,2));
%         end
%         imagesc(data);hc=colorbar;
%         return;
%     end
    [mids,seps,ii,kk]=midkonf2d(N);
%     datums=ones(length(seps),length(mids))*0;%NaN;
    di=diff(unique(N.elec(:,1)));
    del=median(di);
%     ddel=median(diff(unique(mids)));
%     if (ddel<del)&(ddel>=1), del=ddel; end
%     dd=del-floor(del);
%     if (floor(dd*20)~=dd*20)|(length(unique(di))>size(N.elec,1)/5), del=1; end
    cla reset;
    set(gca,'XLim',[min(mids)-del max(mids)+del]);
    set(gca,'YLim',[0.5 length(seps)+0.5]);
    set(gca,'YTickMode','auto','XTickMode','auto')
    if cauto,
        cmin=min(feld(isfinite(feld)));cmax=max(feld(isfinite(feld)));        
    end
    if cmin<0, mm=max(abs([cmin cmax]));
        cmin=-mm;cmax=mm; end
    if clog,
        feld=log10(feld);
        cmin=log10(cmin);
        cmax=log10(cmax);
    end
    if cmin>=cmax, cmin=cmax-0.01*abs(cmax)-1e-5; end
    caxis([cmin cmax]);
    switch cmap
        case 2, colormap(b2r);
        case 3, colormap(hot);
        case 4, colormap(gray);
        case 5, colormap(jet);
        case 6, colormap(cool);
        otherwise, colormap(jet);
    end
    cmap=colormap;lcm=length(cmap);
    for l = 1:length(feld),
        if kk(l)*ii(l)>0, 
            xp=mids(ii(l))+[1 -1 -1 1]*del/2;
            zp=kk(l)+[1 1 -0.95 -0.95]/2;
            if isfinite(feld(l)),
                cind=1+round((feld(l)-cmin)/(cmax-cmin)*lcm);
                if cind>lcm, cind=lcm; end
                if cind<1, cind=1; end
                patch(xp,zp,cmap(cind,:),'EdgeColor',cmap(cind,:),'UserData',l);
            end
        end
    end
    if clog, feld=10.^feld; end
    if isfield(MAL,'equal')&&(MAL.equal>0), axis('equal','tight'); end
    %alpha(1-isnan(datums));
%     cax=caxis;cax(1)=cax(1)-(cax(2)-cax(1))/62;caxis(cax);
%     cmap=colormap;cmap(1,:)=1;colormap(cmap);
    set(gca,'XAxisLocation','top','YDir','reverse','XLim',...
        [min(N.elec(:,1)) max(N.elec(:,1))]);
%     yl=get(gca,'YLim');set(gca,'YTick',yl(1):yl(end)); % Test to see all
    if xdir, set(gca,'XDir','reverse'); end
    if cbar,
        hc=colorbar('horiz');dar=get(hc,'DataAspectRatio');
        set(hc,'DataAspectRatio',dar.*[1 64 1]);
%         old system of alpha-shading        
%         xl=get(hc,'XLim');xl(1)=xl(1)+diff(xl)/64;set(hc,'XLim',xl);
        xt=get(hc,'XTick');
        xtl=num2strcell(rndig(xt));
        if (min(feld(:))>0)&&(clog),
            xtl=num2strcell(round(10.^xt));
            fi=find(xt<1);
            for i=1:length(fi), xtl{fi(i)}=num2str(0.01*round(100*10.^xt(fi(i)))); end
        end
        if isfield(MAL,'canot')&&ischar(MAL.canot),
            xtl{end-1}=MAL.canot;
        end
        set(hc,'XTickMode','manual','XTickLabel',xtl);
    end
    % if min(feld(:))>0,
    %     set(get(hc,'XLabel'),'String','\rho_a in \Omega m');
    % else
    %     set(get(hc,'XLabel'),'String','\Delta\rho_a in %');    
    % end
    ytl='';
    %    0 1 2 3 4 * 100
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
%     xl=get(gca,'XTickLabel');
%     xl(end-1,1:3)='x/m';
    xt=get(gca,'XTick');
    xl={};for iii=1:length(xt), xl{iii}=num2str(rndig(xt(iii))); end
    xl{end-1}='x/m';
    set(gca,'XTickMode','manual','XTickLabel',xl);
%     set(get(gca,'XLabel')
    hold on
    plot(N.elec(:,1),ones(size(N.elec,1),1)*0.5,'k.','MarkerSize',1);
    hold off
    set(gca,'FontSize',8);
end % surface data
set(gca,'TickLength',get(gca,'TickLength')/5);
N=NN;
el0=find(N.elec(:,2)==0);
nn.elec=N.elec(el0,:);
% fsurf=findmess(N,el0,el0);
% number of boreholes
nbh=unique(N.elec(N.elec(:,2)>0,1));
for nb=1:length(nbh),
  el=find((N.elec(:,1)==nbh(nb)).*N.elec(:,2));
  % surface-hole measurements
  fsh=findmess(N,el,el0);
  if ~isempty(fsh),
      x=unique(N.elec(el0,1));
      z=unique(N.elec(el,2));
      data=ones(length(z),length(x))*NaN;
      for i=1:length(fsh),
        [tf,xloc]=ismember(N.elec(N.a(fsh(i)),1),x);
        [tf,zloc]=ismember(N.elec(N.m(fsh(i)),2),z);
        if isnumeric(xloc)&&isnumeric(zloc)&&(xloc>0)&&(zloc>0), 
            iii=fsh(i);
            if (iii>0)&&(iii<=length(feld)), data(zloc,xloc)=feld(iii); end
        end
      end
      if any(data),
          set(figure(nb+10),'NumberTitle','off','Name',['Borehole ' num2str(nb)]);
          clf;
          iconify(nb+10);
          if clog,
              imagesc(x,z,log10(data));
          else
              imagesc(x,z,data);
          end
          if cauto==0,
              if clog,
                  caxis(log10([cmin cmax]));
              else
                  caxis([cmin cmax]);
              end
          end
          alpha(1-isnan(data));
          axis equal tight
          hc=colorbar('horiz');
          if clog,
              if (max(data(:))<10)||((cauto==0)&&(cmax<10)),
                  xt=num2str(0.01*round(100*10.^str2num(get(hc,'XTickLabel'))));
              else
                  xt=num2str(round(10.^str2num(get(hc,'XTickLabel'))));
              end
              set(hc,'XTickLabel',xt);
          end
          title(['Surface hole for borehole x=' num2str(nbh(nb))]);
          %       xlabel('x in m');
          %       ylabel('z in m');
          xl=get(gca,'XTickLabel');
          xl(end-1,1:5)='x [m]';
          set(gca,'XTickLabel',xl);
          xl=get(gca,'YTickLabel');
          xl(end-1,1:5)='z [m]';
          set(gca,'YTickLabel',xl);
          hold on;plot(N.elec(el,1),N.elec(el,2),'k.','MarkerSize',4); hold off
      end
  end
end
if length(nbh)>1, % crosshole measurements
  el1=find((N.elec(:,1)==nbh(1)));%.*N.elec(:,2));
  el2=find((N.elec(:,1)==nbh(2)));%.*N.elec(:,2));
  fsh=findmess(N,el1,el2);
  z1=unique(N.elec(el1,2));
  z2=unique(N.elec(el2,2));
  data=ones(length(z1),length(z2))*NaN;
  for i=1:length(fsh),
      [tf,loc1]=ismember(N.elec(N.a(fsh(i)),2),z1);
      [tf,loc2]=ismember(N.elec(N.m(fsh(i)),2),z2);
      data(loc2,loc1)=field(fsh(i));
  end
  if find(~isnan(data)),
%       figure(nb+2);    
      set(figure(10),'NumberTitle','off','Name','Cross-borehole ','MenuBar','none');
      iconify(10);
      clf;
      if clog,
          data(data<=0)=NaN;
          imagesc(z1,z2,log10(data));
      else
          imagesc(z1,z2,data);
      end
      if cauto==0,
          if clog,
              caxis(log10([cmin cmax]));
          else
              caxis([cmin cmax]);
          end
      end
      alpha(1-isnan(data));
      axis equal tight
      hc=colorbar;
      if clog,
          if max(data(:))>10,
              xt=num2str(round(10.^str2num(get(hc,'YTickLabel'))));
          else
              xt=num2str(0.01*round(100*10.^str2num(get(hc,'YTickLabel'))));
          end
          set(hc,'YTickLabel',xt);
      end
      xlabel('z1 in m');
      ylabel('z2 in m');
      title('Cross-borehole 1-2');
  end
end
if clog, cmin=10.^cmin;cmax=10.^cmax; end
if nargout<1, mids=[cmin cmax]; end

function fi=findmess(N,el1,el2)
el1=[el1(:);0];
el2=[el2(:);0];
a1=ismember(N.a,el1);
fi=find(N.b);a1(fi)=a1(fi)&ismember(N.b(fi),el1);
a2=ismember(N.a,el2);
fi=find(N.b);a2(fi)=a2(fi)&ismember(N.b(fi),el2);
m1=ismember(N.m,el1);
fi=find(N.n);m1(fi)=m1(fi)&ismember(N.n(fi),el1);
m2=ismember(N.m,el2);
fi=find(N.n);m2(fi)=m2(fi)&ismember(N.n(fi),el2);
fi=find(((a1+a2).*(m1+m2)>0)+((a1+m1).*(a2+m2)>0)==2);
