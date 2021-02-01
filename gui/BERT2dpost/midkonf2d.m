function [mids,seps,ii,kk,midpoint,dm]=midkonf2d(N)

% MIDKONF2D - midpoint and konfiguration of data
% [mids,seps,ii,kk]=midkonf2d(N)
% mids - midpoint/reference point of datum
% konf - konfiguration = 
%   type*10000 + dipole_length*100 + separation
% in relation to minimum electrode spacing
% type = 0 - pole-pole
%        1 - pole-dipole forward
%        2 - pole-dipole reverse
%        3 - Wenner (or other wenner)
%        4 - Schlumberger
%        5 - dipole-dipole (or other gamma)
%        6 - dipole-dipole (or other gamma)
%        7 - dipole-dipole (or other gamma)

del=diff(unique(N.elec(:,1)));
%dx=min(del(find(del>1e-4)));
dx=median(del(del>1e-4));
dx=round(dx*10000)/10000;
% 2-Punkt = Pol-Pol
sep=abs(N.elec(N.a,1)-N.elec(N.m,1))/dx; % 0 + sep
midpoint=(N.elec(N.a,1)+N.elec(N.m,1))/2; % 2-Punkt
fn=find(N.n>0);
if ~isempty(fn), % 3-Punkt oder MN<<AB
%     midpoint(fn)=(2*N.elec(N.a(fn),1)+N.elec(N.m(fn),1)+...
%         N.elec(N.n(fn),1))/4;
    mn=abs(N.elec(N.m(fn),1)-N.elec(N.n(fn),1));
    midpoint(fn)=(N.elec(N.m(fn),1)+N.elec(N.n(fn),1))/2;
    sep(fn)=abs(N.elec(N.m(fn),1)-N.elec(N.a(fn),1))/dx+...
        10000+10000*abs(N.elec(N.a(fn),1)>N.elec(N.m(fn),1))+...
        100*(mn/dx-1);
        %abs(N.elec(N.m(fn),1)-N.elec(N.n(fn),1))+...
end % PD=10000(forward),20000(reverse)
fb=find(N.b>0);
if ~isempty(fb), % 4-Punkt, Vorsicht Dipol-Pol!!!
    ab=abs(N.elec(N.a(fb),1)-N.elec(N.b(fb),1));
    mn=abs(N.elec(N.m(fb),1)-N.elec(N.n(fb),1));
    am=abs(N.elec(N.a(fb),1)-N.elec(N.m(fb),1));
    bn=abs(N.elec(N.b(fb),1)-N.elec(N.n(fb),1));
    an=abs(N.elec(N.a(fb),1)-N.elec(N.n(fb),1));
    bm=abs(N.elec(N.b(fb),1)-N.elec(N.m(fb),1));
    ab=round(ab*10000)/10000;mn=round(mn*10000)/10000;
    am=round(am*10000)/10000;bn=round(bn*10000)/10000;
    bm=round(bm*10000)/10000;
    sep(fb)=abs(N.elec(N.m(fb),1)-N.elec(N.b(fb),1))/dx;
    f1=find(ab~=mn); % CPPC, C-PP-C
    f1b=fb(f1);
    if ~isempty(f1),
        midpoint(f1b)=(N.elec(N.m(f1b),1)+N.elec(N.n(f1b),1))/2;
        spac=abs(N.elec(N.n(f1b),1)-N.elec(N.b(f1b),1));
%         sep(f1)=spac/dx+30000+(3*mn(f1)<ab(f1))*10000;
         abmn3=round((3*mn(f1)-ab(f1))*10000)/10000;
         sep(f1b)=spac./dx+(mn(f1)/dx-1)*100.*(abmn3~=0)+30000+(abmn3<0)*10000;
    end % WE=30000, SL=40000
    f3=fb(find((ab>mn).*(abs(am-bn)>1e-4))); % find MUST stay here!
    if ~isempty(f3),  % C---PP--C (gradient type)
%         sep(f3)=70000+mn(f3)/dx+(ab(f3)/dx/2-1)*100; 
%         sep(f3)=70000+mn(f3)/dx*1000+ab(f3)/dx+am(f3)*10+bn(f3)*100;
%         sep(f3)=70000+mn(f3)/dx*1000+ab(f3)/dx+min(am(f3),bn(3))/dx*10-abs(am(f3)-bn(f3))/dx*100; 
%         sep(f1)=70000+(mn(f1)/dx-1)*100+ab(f1)/dx;
%         sep(f1)=70000+(mn(f1)/dx-1)*100+ab(f1)/dx; % assume there are no additional wenner/schlumberger
        mab=(N.elec(N.a(f1),1)+N.elec(N.b(f1),1))/2;
        mmn=(N.elec(N.m(f1),1)+N.elec(N.n(f1),1))/2;        
%         sep(f1)=70000+(ab(f1)/dx+1)*100+mmn-mab;
        sep(f1)=70000+(ab(f1)/dx+1)*100+abs(mmn-mab);
        fi=find(mmn-mab<0);sep(f1(fi))=sep(f1(fi))+10000;
        midpoint(f3)=(N.elec(N.m(f3),1)+N.elec(N.n(f3),1))/2;
%         dx=dx/2;
    end
    f2=find(ab==mn); % CPCP,CCPP, CC-PP
    if ~isempty(f2), % CCPP=50000
        spac=min([am(f2) bn(f2) bm(f2) an(f2)],[],2);
%         sep(fb(f2))=spac/dx+50000+(ab(f2)/dx-1)*10; 
%         sep(fb(f2))=spac./ab(f2)+50000+(ab(f2)/dx-1)*100; 
        sep(fb(f2))=spac./dx+50000+(ab(f2)/dx-1)*100; 
    end
%     midpoint(fb(f2))=(N.elec(N.b(fb(f2)),1)+N.elec(N.m(fb(f2)),1))/2;
    midpoint(fb(f2))=(N.elec(N.b(fb(f2)),1)+N.elec(N.m(fb(f2)),1)+...
        N.elec(N.a(fb(f2)),1)+N.elec(N.n(fb(f2)),1))/4;
end
sep=round(sep*10)/10;
% sep(find(sep==0))=9999; % Circulating=9999
fi=find(fix(sep/10000)==4); % single wenner/schlumberger=circulate!
if length(fi)==1, sep(fi)=70000+mn(fi); end
if max(sep)<0, sep=abs(sep); end
seps=unique(sep);
% [aa,bb]=meshgrid(sep,seps);
% [kk,jj]=find((aa-bb)==0);
[tf,kk]=ismember(sep,seps);
midpoint=round(midpoint*1000)/1000;
mids=unique(midpoint); %alt!
dm=min(diff(mids));
% if dm>5*dx, %possibly 1d
%     dm=dx; end
if isempty(dm)||(dm<0.0499)&&(max(mids)-min(mids)>10), dm=1; end
if (dm<0.001)&&(max(mids)-min(mids)>0.1), 
    dm=0.001;
    mids=round(mids*1000)/1000;
    midpoint=round(midpoint*1000)/1000;
end
% if isfield(N,'eind'), dm=10^fix(log10(dm)); end
dm=round(dm*1000)/1000;
if isempty(dm)||(dm<=0), dm=1; end
mul=dm/2;

% [tf,ii]=ismember(round(midpoint/mul)*mul,round(mids/mul)*mul);
mids=round((min(mids):dm:max(mids))/mul)*mul;
[tf,ii]=ismember(round(midpoint/mul)*mul,mids);
% [aa,bb]=meshgrid(midpoint,mids);
% [ii,jj]=find((aa-bb)==0);

