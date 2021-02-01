function modelexport2d(outfile,M,x,z,Cov,midp,xz,isip)

% MODELEXPORT2D - export 2d-model o ASCII file
% modelexport2d(outfile,M,x,z);
% modelexport2d(outfile,M,x,z,Coverage);
% modelexport2d(outfile,M,x,z,[],1); % save by midpoint

error(nargchk(4,9,nargin));
if nargin<8, isip=0; end
if nargin<7, xz=zeros(size(x)); end
if nargin<6, midp=0; end
if nargin<5, Cov=[]; end
fid=fopen(outfile,'w');
if fid<0, error('File cannot be opened!'); end
xm=(x(1:end-1)+x(2:end))/2;
zm=(z(1:end-1)+z(2:end))/2;

st=sprintf('#%s\t%s\t%s\t%s\t%s','x1/m','x2/m','z1/m','z2/m','rho/Ohmm');
if isequal(size(M),size(Cov)), 
    if isip, 
        st=[st ' phase']; 
    else
        st=[st ' coverage']; 
    end
end
if midp, % xm zm rho
    fprintf(fid,'#%s\t%s\t%s\r\n','xm/m','zm/m','rho/Ohmm');
    for i = 1:min(length(zm),size(M,2)),
        for j = 1:min(length(xm),size(M,1)),
            zz=zm(i);if xz(j)>0, zz=(xz(j)+xz(j+1))/2-zm(i); end
            fprintf(fid,'%.2f\t%.2f\t%.2f\r\n',xm(j),zz,M(j,i));
        end
    end    
else % x1 x2 z1 z2 rho
    fprintf(fid,'%s\r\n',st);
    if isequal(size(M),size(Cov)),
        for k = 1:size(M,2),
            for i = 1:size(M,1),
		if xz(i)==0,
                  fprintf(fid,'%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\r\n',...
                      x(i),x(i+1),z(k),z(k+1),M(i,k),Cov(i,k));
		else
                  fprintf(fid,'%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.3f\r\n',...
                      x(i),x(i+1),xz(i)-z(k),xz(i)-z(k+1),...
                      xz(i+1)-z(k),xz(i+1)-z(k+1),M(i,k),Cov(i,k));
		end
            end
        end    
    else
        for k = 1:size(M,2),
            for i = 1:size(M,1),
		if xz(i)==0,
                  fprintf(fid,'%.2f\t%.2f\t%.2f\t%.2f\t%.2f\r\n',...
                      x(i),x(i+1),z(k),z(k+1),M(i,k));
		else
                  fprintf(fid,'%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\r\n',...
                      x(i),x(i+1),xz(i)-z(k),xz(i)-z(k+1),...
                      xz(i+1)-z(k),xz(i+1)-z(k+1),M(i,k));
		end
            end
        end    
    end
end
fclose(fid); 
