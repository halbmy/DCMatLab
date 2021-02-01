function Data = loadradiccyl(filename,rad)

if nargin<2, rad=0.15; end
fid=fopen(filename,'r');
header=fgetl(fid);
A=mytextscan(fid,'%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
fclose(fid);
%%
cycle=A{9};
unc=unique(cycle);
for i=1:length(unc),
    fi=find(cycle==unc(i));
%     phi=A{2}(fi)*pi/180;AA=[cos(phi)*rad sin(phi)*rad -A{1}(fi)/100];
%     phi=A{4}(fi)*pi/180;BB=[cos(phi)*rad sin(phi)*rad -A{3}(fi)/100];
%     phi=A{6}(fi)*pi/180;MM=[cos(phi)*rad sin(phi)*rad -A{5}(fi)/100];
%     phi=A{8}(fi)*pi/180;NN=[cos(phi)*rad sin(phi)*rad -A{7}(fi)/100];
%     N=abmn2n(AA,BB,MM,NN);
    N=abmn2n([-A{1}(fi)/100 A{2}(fi)],[-A{3}(fi)/100 A{4}(fi)],[-A{5}(fi)/100 A{6}(fi)],[-A{7}(fi)/100 A{8}(fi)]);
    N.elec(:,3)=N.elec(:,1);
    phi=N.elec(:,2)*pi/180;N.elec(:,1)=cos(phi)*rad;N.elec(:,2)=sin(phi)*rad;
    N.rho=A{12}(fi)./A{17}(fi);
    N.ip=A{13}(fi);
    N.err=A{14}(fi)/100;
%     N.time=A{16}(fi);
    if length(unc)>1, Data{i}=N; else Data=N; end
end
