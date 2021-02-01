function modelexport3d(outfile,M,x,y,z,Cov)

% MODELEXPORT3d - Export 3d model to ASCII file
% modelexport3d(modfile,M,x,y,z[,Cov])
% modelexport3d(modfile,Model[,Cov])

if isstruct(M),
    iscov=0;
    if nargin>2, 
        Cov=x;
        iscov=isequal(size(Cov),size(M.M));
    end
    x=M.x;
    y=M.y;
    z=M.z;
    M=M.M;
else
    iscov=(nargin>5)&&isequal(size(Cov),size(M));
end

% ss='%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f';
ss='%g\t%g\t%g\t%g\t%g\t%g\t%.2f';
if iscov, ss=[ss '\t%.2f']; end
ss=[ss '\r\n'];
fid=fopen(outfile,'w'); 
fprintf(fid,'#%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n',...
    'x1','x2','y1','y2','z1','z2','rho_c','(cov)');
for k = 1:size(M,3),
    z1=z(k);z2=z(k+1);
    for j = 1:size(M,2),
        y1=y(j);y2=y(j+1);
        for i = 1:size(M,1),
            x1=x(i);x2=x(i+1);
            if iscov,
                fprintf(fid,ss,x1,x2,y1,y2,z1,z2,M(i,j,k),Cov(i,j,k));
            else
                fprintf(fid,ss,x1,x2,y1,y2,z1,z2,M(i,j,k));
            end
        end
    end
end     
fclose(fid);
