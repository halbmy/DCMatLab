dx=2.5;
xa=min(N.elec(:,1));
xe=max(N.elec(:,1));
ya=min(N.elec(:,2));
ye=max(N.elec(:,2));
[XE,YE]=meshgrid(xa:dx:xe,ya:dx:ye);
ZE=griddata(N.elec(:,1),N.elec(:,2),N.elec(:,3),XE,YE);
clf
surface(XE,YE,ZE);
shading flat
hold on
plot3(N.elec(:,1),N.elec(:,2),N.elec(:,3),'k.')
hold off
colorbar
text(2042.38,1730,'Höhe in m');
xlabel('Rechtswert - 3595000');
ylabel('Hochwert - 5720000');
return
for i=1:length(N.zweid),
    NN=N.zweid{i};
    le=size(NN.elec,1);
    [i le length(find(NN.elec(:,2)==NaN))];
%     fi=find(NN.    
end
