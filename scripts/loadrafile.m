global datfile
[fpath,fname,fext]=fileparts(datfile);
rafile=strrep(datfile,fext,'.ra');
A=load(rafile);
un=unique(A(:,1));
un(un==1)=[];
hold on
for i=1:length(un),
    fi=find(A(:,1)==un(i));
    xx=A(fi,2);
    zz=A(fi,3);
    set(line(xx,zz),'Color','black','LineWidth',1.5);
end
hold off
    