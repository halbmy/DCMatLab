function savecirc2dfile(fname,N)

% SAVECIRC2DFILE - Save INV2D file format (circular)
% savecirc2dfile(filename,N)

fid=fopen(fname,'w');
s=1;
if fid<0,
    error('File not found!');
end
%newline='\r\n';
newline='\n';
fprintf(fid,['%d' newline],size(N.elec,1));
%fprintf(fid,'# Positions(x,z) for all electrodes\r\n');
ss='';
for l=1:size(N.elec,2), ss=[ss '\t%g']; end
ss=[ss newline];
ss(1:2)=[];
fprintf(fid,ss,N.elec');
fprintf(fid,['%d' newline],length(N.r));
%fprintf(fid,'# Electrode numbers A B M N (0=inf), rhoa, error, ...\r\n');
mess=[N.a(:) N.b(:) N.m(:) N.n(:) N.r(:)];
ss=['%3d  %3d  %3d  %3d  %.3f'];
doku='#A B M N rho_a';
if isfield(N,'err'),
    mess=[mess N.err(:)];
    ss=[ss '  %.3f'];
    doku=[doku ' error'];
end
% if isfield(N,'i'),
%     mess(1:length(N.i),end+1)=N.i(:);
%     ss=[ss '  %.4f'];
%     doku=[doku ' I/A'];
% end
% if isfield(N,'u'),
%     mess(1:length(N.u),end+1)=N.u(:);
%     ss=[ss '  %.4f'];
%     doku=[doku ' U/V'];
% end
% if isfield(N,'ip')&&(~isempty(find(N.ip))),
%     mess(1:length(N.ip),end+1)=N.ip(:);
%     ss=[ss '  %.3f'];
%     doku=[doku ' IP'];
% end
ss=[ss newline];
doku=[doku newline];
%fprintf(fid,doku);
fprintf(fid,ss,mess');
fclose(fid);
