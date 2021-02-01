function Data1=extractsounding(Data)

% EXTRACTSOUNDING - Extract median sounding from 2d data set
% Data1 = extractsounding(Data)

am=Data.a-Data.m;
bm=Data.b-Data.m;
[unabm,II,JJ]=unique([am bm],'rows');
medabmn=(Data.a+Data.m+Data.b+Data.n)/4;
medel=floor(size(Data.elec,1)/2);
Data1=[];
Data1.elec=Data.elec;
Data1.r=zeros(length(II),1);
for i=1:length(II),
    fi=find(JJ==i);
    Data1.r(i)=median(Data.r(fi));
%     [mi,fm]=min(abs(Data.a(fi)-median(Data.a(fi))));
%     [mi,fm]=min(abs(medabmn(fi)-medel));
%     fm=fm(1);
    fm=1;
    Data1.a(i)=Data.a(fi(fm));Data1.b(i)=Data.b(fi(fm));
    Data1.m(i)=Data.m(fi(fm));Data1.n(i)=Data.n(fi(fm));
end