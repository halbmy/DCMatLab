function writetom(Shot,filename)

% WRITETOM - Write *.tom file
% writetom(Shot,filename)

if nargin<1, error('must specify Shot!'); end
if nargin<2, filename='test.tom'; end
A=zeros(length(Shot.t),8);
A(:,1)=Shot.t*1000;
ges=0;
for i=1:length(Shot.ns),
    le=length(Shot.nx{i});
    A(ges+1:ges+le,[3 5])=repmat(Shot.pos(Shot.ns{i},:),le,1);
    A(ges+1:ges+le,[6 8])=Shot.pos(Shot.nx{i},:);
    ges=ges+le;
end
fid=fopen(filename,'w');
fprintf(fid,'%.2f\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\r\n',A');
fclose(fid);
return
topo=load('ea\profila.topo');
Shot.pos(:,2)=interp1(topo(:,1),topo(:,2),Shot.pos(:,1),'linear','extrap');
Shot.pos=massbandkorr(Shot.pos);