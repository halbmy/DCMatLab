function createres2dinvfile(typ,fname,ne,x0,dx,n0,nmax)

% CREATERES2DINVFILE - create 2D Loke & Barker Data File
% createres2dinvfile(typ,fname,ne,x0,dx,n0,nmax)
% typ - type of array
%       1 - Wenner (cppc)
%       2 - Pole-Pole (cp)
%       3 - Dipole-Dipole (cc-pp)
%       4 - Wenner-beta (ccpp)
%       5 - Wenner-gamma (cpcp)
%       6 - Pol-Dipol (c-pp)
%       7 - Schlumberger (c-pp-c)
%       8 - Median Gradient (for ccppc)
% fname - filename to write in[stdout]
% ne - number of Electrodes[42]
% x0 - position of first electrode[0]
% dx - electrode spacing[1]
% n0 - minimum separation[1]
% nmax -maximum seperation[8]

if nargin<1, error('Type of measuremente required!'); end
if nargin<7, nmax=8; end % 8
if nargin<6, n0=1; end
if nargin<5, dx=1; end
if nargin<4, x0=0; end
if nargin<3, ne=42; end % 42
if nargin<2, fname=''; end
if (typ>7)|(typ<1), error('Type of measurement unknown'); end

rho=100;
nntyp=[3 3 4 3 3 4 4 10]; %No. of columns
arrayname={'Wenner','Pole-Pole','Dipol-Dipole','Wenner-beta','Wenner-gamma','Pol-Dipol','Schlumberger','Median Gradient'};
out=[];
for n=n0:nmax,
    first=1;
    switch typ,
    case {1,4,5}, % WENNER
      while first+n*3<=ne,
          out=[out;[x0+dx*(first-1) dx*n rho]];
          first=first+1;
      end
    case 2, % POL_POL
      while first+n<=ne,
          out=[out;[x0+dx*(first-1) dx*n rho]];
          first=first+1;
      end
    case 3, % DIPOL-DIPOL (std)
      while first+n+2<=ne,
          out=[out;[x0+dx*(first-1) dx n rho]];
          first=first+1;
      end
    case 6, % POL_DIPOL (std)
      while first+n+1<=ne,
          out=[out;[x0+dx*(first-1) dx n rho]];
          first=first+1;
      end
    case 7,
      while first+2*n+1<=ne,
          out=[out;[x0+dx*(first-1) dx n rho]];
          first=first+1;
      end
    case 8,
    end
end

if isempty(fname),
    fid=1;    
    fin='\n';
else
    fid=fopen(fname,'w');
    fin='\r\n';
end
ss='';
for n=1:nntyp(typ), ss=[ss '%.1f ']; end
ss=[ss fin];
data=size(out,1);

fprintf([ num2str(data) ' ' arrayname{typ} ' ' num2str(ne) ' Electrodes, x0=' num2str(x0) ' dx=' num2str(dx) ' n=' num2str(n0) '-' num2str(nmax) fin]);
fprintf(fid,[ arrayname{typ} ' ' num2str(ne) ' Electrodes, x0=' num2str(x0) ' dx=' num2str(dx) ' n=' num2str(n0) '-' num2str(nmax) fin]);
fprintf(fid,['0' fin]); % xgrid
fprintf(fid,[num2str(typ) fin]); % type
fprintf(fid,[num2str(data) fin]); % no. of data
fprintf(fid,['0' fin]); % x0=first electrode
fprintf(fid,['0' fin]); % IP present
fprintf(fid,ss,out'); %DATA
fprintf(fid,['0' fin '0' fin '0' fin '0' fin '0' fin]); % END

if fid~=1, fclose(fid); end
