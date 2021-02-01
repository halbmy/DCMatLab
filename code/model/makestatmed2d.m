function B=makestatmed(Mod,range,nreal,randstart)
%GSLIB2DC3D
%wandelt mit lusim (GSLIB) generierte zufallsverteilte Medien (mod.
%Parameter in Spaltenform) um in dc3dinf (T.Guenther) Format (*.mod)
%fuer jede realisation wird neues file erstellt
%ji 01/2006

%Eingabe der statistischen Parameter fuer Verteilung
%Mittelwert und Standardabweichung von lg(Leitfähigkeit)
%Umrechnung in rho erfolgt ggf unten

% nx=101;ny=101;nz=10;

if nargin<4, randstart=round(rand*5000)*2+1; end
if nargin<3, nreal=1; end
if nargin<2, range=5; end
while length(range)<3, range(end+1)=range(end); end

dx=median(diff(Mod.x));
nx=floor((max(Mod.x)-min(Mod.x))/dx)+1;

if isfield(Mod,'y'), % 3d model
    dy=median(diff(Mod.y));
    ny=floor((max(Mod.y)-min(Mod.y))/dy)+1;
else
    dy=dx;ny=1;
end
dz=min(diff(Mod.z));
nz=floor((max(Mod.z)-min(Mod.z))/dz)+1;

%SGSIM parameter file (*.par) einlesen
filename='sgsim.par';%tempname;
outfile='sgsim.out';
fid=fopen(filename,'w');
fprintf(fid,'                  Parameters for SGSIM\n');
fprintf(fid,'                  ********************\n');
fprintf(fid,'\n');
fprintf(fid,'START OF PARAMETERS:\n');
fprintf(fid,'\n');
fprintf(fid,'1  2  0  3  5  0              \\  columns for X,Y,Z,vr,wt,sec.var.\n');
fprintf(fid,'-1.0       1.0e21             \\  trimming limits\n');
fprintf(fid,'0                            \\ transform the data (0=no, 1=yes)\n');
fprintf(fid,'sgsim.trn                     \\  file for output trans table\n');
fprintf(fid,'0                             \\  consider ref. dist (0=no, 1=yes)\n');
fprintf(fid,'histsmth.out                  \\  file with ref. dist distribution\n');
fprintf(fid,'1  2                          \\  columns for vr and wt\n');
fprintf(fid,'0.0    15.0                   \\  zmin,zmax(tail extrapolation)\n');
fprintf(fid,'1       0.0                   \\  lower tail option, parameter\n');
fprintf(fid,'1      15.0                   \\  upper tail option, parameter\n');
fprintf(fid,'1                             \\debugging level: 0,1,2,3\n');
fprintf(fid,'sgsim.dbg                     \\file for debugging output\n');
fprintf(fid,'%s                     \\file for simulation output\n',outfile);
fprintf(fid,'%d                         \\number of realizations to generate\n',nreal);
fprintf(fid,'%d   0.0    %.2f             \\nx,xmn,xsiz\n',nx,dx);
fprintf(fid,'%d    0.0    %.2f              \\ny,ymn,ysiz\n',ny,dy);
fprintf(fid,'%d    0    %.2f                \\nz,zmn,zsiz\n',nz,dz);
fprintf(fid,'%d                         \\random number seed\n',randstart);
fprintf(fid,'0     8                       \\min and max original data for sim\n');
fprintf(fid,'12                            \\number of simulated nodes to use\n');
fprintf(fid,'1                             \\assign data to nodes (0=no, 1=yes)\n');
fprintf(fid,'1     3                       \\multiple grid search (0=no, 1=yes),num\n');
fprintf(fid,'0                             \\maximum data per octant (0=not used)\n');
fprintf(fid,'10.0  10.0  10.0              \\maximum search radii (hmax,hmin,vert)\n');
fprintf(fid,'0.0   0.0   0.0              \\angles for search ellipsoid\n');
fprintf(fid,'0     0.60   1.0              \\ktype: 0=SK,1=OK,2=LVM,3=EXDR,4=COLC\n');
fprintf(fid,'../data/ydata.dat             \\  file with LVM, EXDR, or COLC variable\n');
fprintf(fid,'4                             \\  column for secondary variable\n');
fprintf(fid,'1    0                        \\nst, nugget effect\n');
fprintf(fid,'1    1.0  90.0   0.0   0.0     \\it,cc,ang1,ang2,ang3\n');
fprintf(fid,'%.2f  %.2f  %.2f              \\a_hmax, a_hmin, a_vert\n',range(1),range(2),range(3));
fclose(fid);
dos('SGSIM2d.exe');

data=textread(outfile,'','headerlines',3);
% data=(10.^(data.*sigma+mue)).^(-1);%falls lg simga simuliert wurde (log-norm-vtlg)

nsamples=floor(nx*ny*nz); %anzahl der datenpunkte pro realisation
nreal=floor(length(data)/nsamples); %anzahl der realisationen
a=ceil(sqrt(ceil(sqrt(nreal)/0.5))); b=ceil(nreal/a);%anzahl der zeilen und spalten fuer subplots

if nreal>1,
    A1=squeeze(zeros(nx,ny,nz));
    for l=1:nreal,
        dat=data((l-1)*nsamples+1:l*nsamples);
        A1(:)=dat;
        A{l}=A1;
    end
else
    A=squeeze(zeros(nx,ny,nz));
    A(:)=data(1:nsamples);
end
% delete(filename);
%%
x1=Mod.x(1)+dx/2+(0:nx-1)*dx;
z1=Mod.z(1)+dz/2+(0:nz-1)*dz;
[Z1,X1]=meshgrid(z1,x1);
xx=Mod.x(1:end-1)+diff(Mod.x)/2;
zz=Mod.z(1:end-1)+diff(Mod.z)/2;
if nreal>1,
    for l=1:nreal,
       B{l}=interp2(Z1,X1,A{l},zz(:)',xx(:),'linear');
    end
else
    B=interp2(Z1,X1,A,zz(:)',xx(:),'linear');
end
% imagesc(Mod.x,Mod.z,B);colorbar;axis equal tight
