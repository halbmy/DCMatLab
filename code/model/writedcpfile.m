function writedcpfile(N,filename,comment)

% WRITEDCPFILE - write data structure to DCP file for em1dinv use
%                uses source type 24 (arbitrary xyz electrode locations)
% writedcpfile(Data,filename[comment])
% Data consists of fields: elec (electrode positions

if nargin<1, error('Specify file name and data structure'); end
if nargin<2, filename='sond.dcp'; end
if nargin<3, comment='Sounding'; end
if isstruct(N), % normal data structure
    A=zeros(length(N.a),14);
    A(:,1)=N.elec(N.a,1); % x coordinate
    A(:,4)=N.elec(N.m,1);
    A(:,7)=N.elec(N.n,1);
    A(:,10)=N.elec(N.b,1);
    if size(N.elec,2)>1, % depth
        A(:,3)=N.elec(N.a,end);
        A(:,6)=N.elec(N.m,end);
        A(:,9)=N.elec(N.n,end);
        A(:,12)=N.elec(N.b,end);
    end
    if size(N.elec,2)==3, % also y given
        A(:,2)=N.elec(N.a,2);
        A(:,5)=N.elec(N.m,2);
        A(:,8)=N.elec(N.n,2);
        A(:,11)=N.elec(N.b,2);
    end
    A(:,13)=N.r;
    A(:,14)=0.05;
else % ab2/mn2/data array
    if size(N,2)==3, %only ab2/mn2/data->type 22
        fid=fopen(filename,'w');
        fprintf(fid,'%s\r\n',comment);
        fprintf(fid,'%d\r\n',22);
        fprintf(fid,'%d\t%d\t%d\t%d\r\n',1,1,1,1);
        for i=1:size(N,1),
           fprintf(fid,'%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.2f\n',...
               -N(i,1),-N(i,2),N(i,2),N(i,1),N(i,3),0.01);
        end
        fclose(fid);
        return;
    else % all positions
        
    end
end
fid=fopen(filename,'w');
fprintf(fid,'%s\r\n',comment);
fprintf(fid,'%d\r\n',24);
fprintf(fid,'%d\t%d\t%d\t%d\r\n',1,1,1,1);
formstr='%.2f';
for i=2:size(A,2), formstr=[formstr '\t%.2f']; end
formstr=[formstr '\r\n'];
fprintf(fid,formstr,A');
fclose(fid);