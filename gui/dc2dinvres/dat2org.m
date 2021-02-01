% du=N.a;N.a=N.m;N.m=du;du=N.b;N.b=N.n;N.n=du;
all=[N.a(:) N.b(:) N.m(:) N.n(:)];
allsort=sortrows(all);
olda=0;oldb=0;
l=0;n=0;org=zeros(round(length(N.a)/4),10);
for i=1:length(N.a),
    n=n+1;
    if (allsort(i,1)~=olda)||(allsort(i,2)~=oldb)||(n>4),
        n=1;
        l=l+1;
    end
    org(l,1)=allsort(i,1);
    org(l,2)=allsort(i,2);
    org(l,n*2+1)=allsort(i,3);
    org(l,n*2+2)=allsort(i,4);
    olda=allsort(i,1);oldb=allsort(i,2);
end
% fi=find(org(:,5)==0);org(fi,:)=[];
fprintf('Effectivity %d%%\n',round(length(N.a)/size(org,1)/4*100));
fid=fopen('test4.org','w');
fprintf(fid,'10\r\nSCHLUMBERGERINVERS\r\n');
fprintf(fid,'%2d\t%2d\t%2d\t%2d\t%2d\t%2d\t%2d\t%2d\t%2d\t%2d\r\n',org');
fclose(fid);