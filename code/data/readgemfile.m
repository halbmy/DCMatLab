function N=readgemfile(gemfile)

% READGEMFILE - read sensinv3d file
% [fname,fpath]=uigetfile('*.gem');
% gemfile=fullfile(fpath,fname);

fid=fopen(gemfile,'r');
if fid<0, error('GEM-File does not exist!'); end
A=fscanf(fid,'%d %f %f %f',[4 Inf])';
fclose(fid);
N.elec=A(:,2:4);
senfile=strrep(gemfile,'.gem','.sen');
fid=fopen(senfile,'r');
if fid<0, error('SEN-File does not exist!'); end
A=fscanf(fid,'%d %d %d %d %d',[5 Inf])'; %Nummer El. beachten!
fclose(fid);
N.a=abs(A(:,2));
N.b=abs(A(:,3));
N.m=abs(A(:,4));
N.n=abs(A(:,5));
uifile=strrep(gemfile,'.gem','.ui');
fid=fopen(uifile,'r');
if fid<0, error('UI-File does not exist!'); end
A=fscanf(fid,'%f %f',[2 Inf])';
fclose(fid);
N.u=A(:,1)/1000;
N.i=A(:,2)/1000;
N.rho=N.u./N.i;
N.k=getkonf(N);
N.r=N.rho.*N.k;