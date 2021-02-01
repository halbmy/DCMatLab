function newelec=massbandkorr(elec)

newelec=elec;
[so,si]=sort(elec(:,1));
for i=2:size(elec,1),
    ii=si(i);ii1=si(i-1);
    dx=elec(ii,1)-elec(ii1,1);
    dz=max(abs(elec(ii,2:end)-elec(ii1,2:end)),[],2);
    newelec(ii,1)=newelec(ii1,1)+sqrt(dx*dx-dz*dz);
end