function N=readfeinv3dfile(fname)

% READINV3DFILE - Read data file in inv3d format
% N = readinv3dfile(filename)
% Format:
% number_of_Electrodes
% x_el1 (y_el1) z_el1
% ...
% x_eln (y_eln) z_eln
% number_of_datapoints
% A_1 B_1 M_1 N_1 R_1 K_1 Err_1
% ...
% A_n B_n M_n N_n R_n K_1 Err_n
% (ABMN = Electrode numbers, 0=infinity)
% R = apparent resistivity
% K = configuration factor
% Err = (relative) Error 

N.elec=[];
fid=fopen(fname,'r');
if fid<0, error('File not found!'); end
ne=fscanf(fid,'%d\n',1);
for n=1:ne,
    el=str2num(fgetl(fid));
    N.elec(n,1:length(el))=el;
end
if size(N.elec,2)<2, N.elec(:,2)=0; end
if size(N.elec,2)<3, N.elec(:,3)=0; end
nm=fscanf(fid,'%d\n',1);
N.a=zeros(nm,1);N.b=N.a;N.m=N.a;N.n=N.a;N.r=N.a;N.k=N.a;N.err=N.a;
for n=1:nm,
    zeile=fgetl(fid);
    mess=str2num(zeile);
    if length(mess)<5, break; end
    N.a(n)=mess(1);N.b(n)=mess(2);N.m(n)=mess(3);N.n(n)=mess(4);
    if length(mess)>4, N.r(n)=mess(5); end
    if length(mess)>5, N.k(n)=mess(6); end
    if length(mess)>6, N.err(n)=mess(7); end
end
fclose(fid);

message(sprintf('%s: %d Measurements with %d Electrodes',fname,nm,ne));