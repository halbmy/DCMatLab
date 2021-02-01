%load ed2dres
mal=struct('canot','\rho in \Omegam ','cmin',25,'cmax',1010);
tripatchmod(Mesh,Mesh.model,Mesh.alfa,mal);
yl=ylim;yl(2)=2;ylim(yl);
%%
stat=(0:5)*40+32;
cols={'black','black','blue','black','black','red'};
for i=1:length(stat),
    t=text(stat(i),12,['S' num2str(i)]);
    set(t,'Color',cols{i},'HorizontalAlignment','center');
end
%%
lw=5;
line(stat(3)+[-1 1]*40,2*[1 1],'Color','blue','LineWidth',lw);
line(stat(6)+[-1 1]*40,2*[1 1],'Color','red','LineWidth',lw);
