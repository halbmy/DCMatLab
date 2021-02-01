function savepro(filename,N)

% SAVEPRO - Save pro-file data
% savepro(filename,N)

if nargin<2, error('2 input arguments required!'); end
if ~isfield(N,'zweid'), error('No profile data in file!'); end
fid=fopen(filename,'w');
[pn,fn]=fileparts(filename);
olddir=pwd;
if ~isempty(pn), cd(pn); end
for n=1:length(N.zweid),
    nn=N.zweid{n};
    if length(nn.a)>0,
        nn.r=N.r(N.nr{n});
        if isfield(N,'names')&&(length(N.names)>=n),
            name=N.names{n};
            if strfind(name,'=')&&isempty(strfind(name,'.dat')), name=[name '.dat']; end
            name=strrep(strrep(name,' ',''),'=','_');
        else
            nu=num2str(n);if length(nu)<2, nu=['0' nu]; end
            name=['file' num2str(n) '.dat'];    
        end
        saveinv2dfile(name,nn);
        xx1=nn.elec(nn.a(1),1);xx2=nn.elec(nn.m(1),1);
        aa=find(N.nr{n}==nn.a(1));
        mm=find(N.nr{n}==nn.m(1));
        %     xz1=N.elec(N.a(aa),1:2);
        fprintf(fid,'%s ',name);
        if isfield(N,'points')&&(length(N.points)>=n),
            fprintf(fid,' %.1f',N.points{n});
        else
            xz1=N.elec(N.a(N.nr{n}(1)),1:2);
            xz2=N.elec(N.m(N.nr{n}(1)),1:2);
            xz0=xz1-(xz2-xz1)*xx1/(xx2-xx1);
            fprintf(fid,name);
            fprintf(fid,' %.1f',[xz0 (xz1+xz2)/2]);
        end
        fprintf(fid,'\n\r');
    end
end
fclose(fid);
chdir(olddir);

function save2dnew(fname,N)
fid=fopen(fname,'w');
if fid<0,
    error('File not found!');
end
fprintf(fid,'%d # Number of electrodes\n',size(N.elec,1));
fprintf(fid,'# Positions(x,z) for all electrodes\n');
ss='';
for l=1:size(N.elec,2), ss=[ss '\t%g']; end
ss=[ss '\n'];
ss(1:2)=[];
fprintf(fid,ss,N.elec');
% for i=1:size(N.elec,1),
%     fprintf(fid,'%s\r\n',num2str(N.elec(i,:)));
% end
fprintf(fid,'%d # Number of measurements\n',length(N.r));
fprintf(fid,'# Electrode numbers A B M N (0=inf), rhoa, error, ...\n');
mess=[N.a(:) N.b(:) N.m(:) N.n(:) N.r(:)];
ss='%3d %3d %3d %3d  %.2f';
if isfield(N,'err'),
    mess=[mess N.err(:)*100];
    ss=[ss '  %.2f'];
end
ss=[ss '\n'];
fprintf(fid,ss,mess');
fclose(fid);
