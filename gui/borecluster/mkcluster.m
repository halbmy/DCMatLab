% Depth   Ca   C  Dichte  F  Gamma H   IP   SP  Fe   K  Kal1   Kal2   Oxy  POR RDeep   RShallow   Si  Susz  Th U  
% A=textread('Staßfurt-B.asc','','headerlines',1);
% spalt=[2 3 4 6 7 10 11 12 14 15 16 18 19 20 21];
% spalt=[2 3 4 6 16 19 21];
% Depth Ca C Dev Dichte  DAz F   Gamma    H IP   SP  Fe  K Kal    Oxy POR    RDeep  RShallow Sal Si Susz  Temp   Th  U 
% A=textread('staßfurt-g1.asc','','headerlines',1);
% spalt=[2 3 5 8 9 12 13 14 15 16 17 20 21 23 24];

[fname,pname]=uigetfile('*.asc;*.ASC','Read ASCII Data File');
if ~isstr(fname), return; end
filename=fullfile(pname,fname);
A=textread(filename,'','headerlines',1);
spalt=1:size(A,2);
sptext=strrep(strrep(num2str(spalt),'   ',' '),'  ',' ');
newtext=inputdlg('Specify columns','Columns',1,{sptext});
spalt=str2num(newtext{1});

fid=fopen(filename,'r');
zeile=fgetl(fid);
fclose(fid);
i=1;col={};
while ~isempty(zeile), 
    [coli,zeile]=strtok(zeile);
    if isempty(coli), break; end
    col{i}=coli;
    i=i+1; 
end

% fi=find((A(:,1)>105)&(A(:,1)<198));A=A(fi,:);
A(A<-990)=NaN;
B=A(:,spalt);
dep=A(:,1);
for i=1:size(B,2),
    aa=B(:,i);ma=max(aa);mi=min(aa);
    fi0=find(isnan(aa));
    fi1=find(~isnan(aa));
    aa(fi0)=interp1(dep(fi1),aa(fi1),dep(fi0),'linear','extrap');
    aa=(aa-mi)/(ma-mi);
    B(:,i)=aa;
end
dep(end+1)=2*dep(end)-dep(end-1);
% [ii,jj]=find(isnan(B));B(ii,:)=[];dep(ii)=[];
% result=fcmcluster(B,struct('c',4));
DI=pdist(B,'seuclidean');
LI=linkage(DI,'ward');
colors={'r','g','b','c','m','y','k'};
for cl=2:7,
    CL=cluster(LI,'maxclust',cl);
    %%
    clf;
    for i=1:length(CL),
        col=colors{CL(i)};if CL(i)==7, col=[1 1 1]/2; end
        patch([-1 1 1 -1]*5,[dep(i) dep(i) dep(i+1) dep(i+1)],col,'EdgeColor',col);
    end
    axis equal tight
    box on
    set(gca,'YDir','reverse','XTick',[]);
    yl=get(gca,'Ylim');dy=50;
    for i=round(yl(1)/dy)*dy:10:round(yl(2)/dy)*dy, line([-5 -3],[i i],'Color','black'); end
    pause;
%     epsprint(gcf,['g1-cl' num2str(cl)],1);
end
