function saveres2dinvfile(outfile,N,field,ipfield)

% SAVERES2DINVFILE -  SAVE Data in 2D Loke&Barker Format
% saveres2dinvfile(outfile,N,N.elec<,field>)
% saves field, otherwise N.r

if nargin<2, error('Too less input arguments!'); end
if nargin<3, field=N.r; end
if nargin<4, ipfield=[]; end

isip=isequal(size(field),size(ipfield));
fid=fopen(outfile,'w');
fprintf(fid,'Mixed Array\r\n');
del=min(diff(N.elec(:,1)));
fprintf(fid,'%.1f\r\n',del);
fprintf(fid,'11\r\n');
fprintf(fid,'0\r\n');
fprintf(fid,'Type of measurement (0=app. resistivity,1=resistance)\r\n');
fprintf(fid,'0\r\n');
data=length(N.r);
fprintf(fid,'%d\r\n',data);
fprintf(fid,'1\r\n');
if isip, fprintf(fid,'1\r\nPhase angle\r\ndeg\r\n0.0,5.0\r\n');
else fprintf(fid,'0\r\n'); end
for l = 1:data,
    n=2;
    s=sprintf('%.1f ',N.elec(N.a(l),:));
    if N.b(l)>0,
        s=[s sprintf('%.1f ',N.elec(N.b(l),:))];
        n=n+1;
    end        
    s=[s sprintf('%.1f ',N.elec(N.m(l),:))];
    if N.n(l)>0,
        s=[s sprintf('%.1f ',N.elec(N.n(l),:))];
        n=n+1;
    end
    if isip, fprintf(fid,'%d %s %.2f %.2f\r\n',n,s,field(l),ipfield(l));
    else fprintf(fid,'%d %s %.2f\r\n',n,s,field(l)); end
end
for l=1:5, fprintf(fid,'0\r\n'); end
fclose(fid);
