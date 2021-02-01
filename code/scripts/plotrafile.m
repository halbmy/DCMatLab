if ~exist('rafile')||isempty(rafile),
  [fpath,fname,fext]=fileparts(datafile);
  rafile=['..\' strrep(datafile,fext,'.ra')];
end
A=load(rafile);
un=unique(A(:,1));
un(un==1)=[];
hold on
for i=1:length(un),
    fi=find(A(:,1)==un(i));
    xx=A(fi,2);
    zz=A(fi,3);
    set(line(xx,-zz),'Color','black','LineWidth',1.5);
end
hold off

