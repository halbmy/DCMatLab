meshname='mesh/bbtest';
% [ELE,NODE,FACE]=loadmesh(meshname);
st=[50 56 7 27 2 3 500;51 54 7 27 1.5 2 500];
rf=[1 2 3 1];
clf
xxx=[43 45 50.6 52.6 63.6 65.6];
zzz=[2.3 2.3 0    0    2.9 2.9];
yyy=[-2 30];
hold on
fill3(xxx([1 end end 1]),yyy([1 1 2 2]),3*[1 1 1 1],'y');
for i=1:length(xxx)-1,
    h=fill3([xxx(i) xxx(i+1) xxx(i+1) xxx(i)],[yyy(1) yyy(1) yyy(2) yyy(2)],...
        [zzz(i) zzz(i+1) zzz(i+1) zzz(i)],'g');
%     set(h,'FaceAlpha',0.5);
end
an1=(xxx(3)-xxx(2))/(zzz(3)-zzz(2));
an2=(xxx(end-2)-xxx(end-1))/(zzz(end-2)-zzz(end-1));
xxx=xxx(2:end-1);zzz=zzz(2:end-1);
zzz(2:end-1)=zzz(2:end-1)+0.5;
xxx(1)=xxx(1)-0.5*an1;xxx(end)=xxx(end)-0.5*an2;
for i=1:length(xxx)-1,
    h=fill3([xxx(i) xxx(i+1) xxx(i+1) xxx(i)],[yyy(1) yyy(1) yyy(2) yyy(2)],[zzz(i) zzz(i+1) zzz(i+1) zzz(i)],'g');
%     set(h,'FaceAlpha',0.5);
end
% for e=1:size(FACE,1),
%     pts=NODE(FACE(e,:),1:3);
%     plot3(pts(rf,1),pts(rf,2),pts(rf,3),'y','MarkerSize',1);
% end
for p=0:6,
    beg=41+p*31;
    ent=beg+30;
%     line(NODE([beg ent],1),NODE([beg ent],2),NODE([beg ent],3));
    for i=0:30,
        plot3(NODE(beg+i,1),NODE(beg+i,2),NODE(beg+i,3)-0.1,'vk','MarkerSize',3);    
    end
end
hold off
[x,y,z]=en2mid(ELE,NODE);
xx=x;yy=y;zz=z;
xlim([min(xx) max(xx)]);
ylim([min(yy) max(yy)]);
zlim([min(zz) max(zz)]);
set(gca,'ZDir','reverse');
view(170,10)
patchacht([50 56 56 50 51 53 53 51],[7 7 27 27 7 7 27 27],[3 3 3 3 1.8 1.8 1.8 1.8]);
% for i=1:size(st,1),
%     patchcube(st(i,1),st(i,2),st(i,3),st(i,4),st(i,5),st(i,6));        
% end
daspect([1 1 1/2]);
zlim([0 4]);
xlabel('x in m');
ylabel('y in m');
zlabel('z in m');
text(50,30,3.5,'Untergrund')
text(50,30,2.5,'Deichkörper')
set(text(55.2,30,2.3,'alter Deich'),'Color',[1 1 1])
text(54.2,30,0.25,'Deckschicht');
% text(50,0,3.5,'Untergrund')
% text(55,2,2.8,'Deichkörper')
% set(text(50,5,2.5,'alter Deich'),'Color',[1 1 1])
% text(50.5,0,0.4,'Deckschicht');