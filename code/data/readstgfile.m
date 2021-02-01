function N = readstgfile(fname)

% READSTGFILE - Read sting device file (*.stg)
% N = readstgfile(filename)
% Format: (# character can be used to comment)
% number_of_Electrodes
% # x y h d # token string for meaning of columns
% x_el1 (y_el1) z_el1
% ...
% x_eln (y_eln) z_eln
% number_of_datapoints
% # a b m n u i # possible token string for meaning of columns
% A_1 B_1 M_1 N_1 R_1 (Err_1)
% ...
% A_n B_n M_n N_n R_n (Err_n)
% (= Electrode numbers, 0=infinity)

fid=fopen(fname,'r');   
if fid<0, error('File not found!'); end
dim=2;
zeile=fgetl(fid);tt=strfind(zeile,'Type');
if strfind(lower(zeile(tt+5:end)),'3d'), dim=3; end
zeile=fgetl(fid);tt=strfind(zeile,'Records');
ndata=str2double(zeile(tt+8:end));
Ax=zeros(ndata,1);Bx=Ax;Mx=Ax;Nx=Ax;N.r=Ax;N.i=Ax;%N.u=Ax;%N.err=Ax;
N.rho=Ax;
zeile=fgetl(fid);
for i=1:ndata,
    zeile=fgetl(fid);
%                      NR ,USER,DATE,HH:MM:SS  ,  R,err,I,rhoa,Ax/y/z,Bx/y/z
%     aa=sscanf(zeile,'%*d,%*s ,%*d,%*d:%*d:%*d, %f,%f,%f,%f,%*s , %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f');
%     if length(aa)<14, N=[];fclose(fid);return; end
%     Ax(i)=aa(5);Bx(i)=aa(8);Mx(i)=aa(11);Nx(i)=aa(14);
%     N.r(i)=aa(4);N.i(i)=aa(3)/1000;N.err(i)=aa(2)/100;N.rho(i)=aa(1);
    ko=strfind(zeile,',');
    if length(ko)<15, N=[];fclose(fid);return; end
    zz=strrep(zeile([ko(4)+1:ko(8)-1 ko(9)+1:ko(20)]),',',' ');
    aa=str2num(zz);
    if length(aa)<14, N=[];fclose(fid);return; end
    Ax(i)=aa(5);Bx(i)=aa(8);Mx(i)=aa(11);Nx(i)=aa(14);
    N.r(i)=aa(4);N.i(i)=aa(3)/1000;
%     N.err(i)=aa(2)/1000;
    N.rho(i)=aa(1);    
end
fclose(fid);
N.elec=unique([Ax;Bx;Mx;Nx]);
[tf,N.a]=ismember(Ax,N.elec);
[tf,N.b]=ismember(Bx,N.elec);
[tf,N.m]=ismember(Mx,N.elec);
[tf,N.n]=ismember(Nx,N.elec);
N.elec(:,2)=0;
N.k=getkonf(N);

[pp,ff,ee]=fileparts(fname);
trnfile=strrep(fname,ee,'.trn');
if exist(trnfile,'file'),
    fid=fopen(trnfile);
    for i=1:3, zeile=fgetl(fid); end
    A=mytextscan(fid,'%f%f');
    fclose(fid);
    N.topo=[A{1} A{2}];
end
