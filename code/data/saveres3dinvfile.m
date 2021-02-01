function saveres3dinvfile(outfile,N,field,raw)

% SAVERES3DINVFILE -  Save Data in RES3DINV Format
%                     (using general electrode array)
% saveres3dinvfile(filename,N[,field[,rawdata]])
% saves field, otherwise N.r

if nargin<2, error('Too less input arguments!'); end
if nargin<3, field=N.r; end
if nargin<4, raw=0; end

data=length(N.r);

fid=fopen(outfile,'w');

if ~raw,
    fprintf(fid,'Mixed Array\r\n');
    fprintf(fid,'1\r\n');
    fprintf(fid,'1\r\n');
    del=min(diff(N.elec(:,1)));
    fprintf(fid,'%.2f\r\n',del);
    del=min(diff(N.elec(:,2)));
    fprintf(fid,'%.2f\r\n',del);
    fprintf(fid,'11\r\n');
    fprintf(fid,'%d\r\n',data);
end
for l = 1:data,
    s=sprintf('%.2f ',N.elec(N.a(l),1:2));
    ne=2;
    if N.b(l)>0,
        s=[s sprintf('%.2f ',N.elec(N.b(l),1:2))];
        ne=ne+1;
    end        
    s=[s sprintf('%.2f ',N.elec(N.m(l),1:2))];
    if N.n(l)>0,
        s=[s sprintf('%.2f ',N.elec(N.n(l),1:2))];
        ne=ne+1;
    end
    fprintf(fid,'%d %s %.3f\r\n',ne,s,field(l));
end
if ~raw, 
    for l=1:5, fprintf(fid,'0\r\n'); end 
end

fclose(fid);
