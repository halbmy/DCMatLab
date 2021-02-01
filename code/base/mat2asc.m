function mat2asc(A,filename)

% MAT2ASC - Save matrix into ascii file
% mat2asc(A,filename)
% mat2asc(filename,A)

if nargin<2, error('Specify matrix and filename'); end
if ischar(A)&&isnumeric(filename),
    du=A;
    A=filename;
    filename=du;
end
fstr='%g';
for i=2:size(A,2), fstr=[fstr '\t%g']; end
fstr=[fstr '\n'];
fid=fopen(filename,'w');
fprintf(fid,fstr,A');
fclose(fid);