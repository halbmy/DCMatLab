function saveinv3dfile(fname,N)

% SAVEINV3DFILE - Save data in inv3d format
% saveinv3dfile(filename,N)
% Format reads as follows:
% number_of_electrodes
%   x_el_1 y_el_1 [z_el_1]
%   ...
%   x_el_E y_el_E [z_el_E]
%   number_of_data
%   a_1 b_1 m_1 n_1 rho_a_1 [Err_1 I_1 U_1 IP_1]
%   ...
%   a_N b_N m_N n_N rho_a_N [Err_1 I_1 U_1 IP_1]
% where x/y/z are the electrode positions and
% a/b/m/n are the used electrode numbers (0=infinity)

if nargin<2, error('Two arguments (filename,datastruct) must be specified!'); end
newline='\r\n';
fid=fopen(fname,'w');
if fid<0, error('Could not open file for writing!'); end
fprintf(fid,'%d # number of electrodes',size(N.elec,1));
fprintf(fid,newline);
fprintf(fid,['# x y z' newline]);
ss='%g';
for l=2:size(N.elec,2), ss=[ss '\t%g']; end
ss=[ss newline];
fprintf(fid,ss,N.elec');
fprintf(fid,'%d # number of data',length(N.a));
fprintf(fid,newline);
mess=double([N.a(:) N.b(:) N.m(:) N.n(:)]);
elst='%4d';sep='\t';%sep=' ';
ss=[elst sep elst sep elst sep elst];%sep='  ';
doku='#a b m n';
if isfield(N,'r')&&(length(N.r)==length(N.a)),
    mess=[mess N.r(:)];
    if min(abs(N.r))>1, ss=[ss sep '%.2f']; else ss=[ss sep '%g']; end
    doku=[doku ' rhoa'];
end
if isfield(N,'rho')&&(length(N.rho)==length(N.a)),
    mess=[mess N.rho(:)];
    if min(abs(N.rho))>1, ss=[ss sep '%.3f']; else ss=[ss sep '%g']; end
    doku=[doku ' R'];
end
if isfield(N,'err')&&(length(N.err)==length(N.a)),
    mess=[mess N.err(:)];
    ss=[ss sep '%g'];%.3f
    doku=[doku ' err'];
end
if isfield(N,'ip')&&(length(N.ip)==length(N.a)),
    mess=[mess N.ip(:)];
    ss=[ss sep '%g'];%.3f
    doku=[doku ' ip'];
end
if isfield(N,'i')&&(length(N.i)==length(N.a)),
    mess=[mess N.i(:)];
    ss=[ss sep '%g'];%.4f
    doku=[doku ' I'];
end
if isfield(N,'u')&&(length(N.u)==length(N.a)),
    mess=[mess N.u(:)];
    ss=[ss sep '%g'];%.4f
    doku=[doku ' U'];
end
if isfield(N,'konf')&&(length(N.konf)==length(N.a)),
    mess=[mess N.konf(:)];
    ss=[ss sep '%g'];%.4f
    doku=[doku ' k'];
end
doku=[doku newline];
fprintf(fid,doku);
ss=[ss newline];
fprintf(fid,ss,mess');
fclose(fid);
