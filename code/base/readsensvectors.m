function S=readsensvectors(Mesh,N)

% READSENSVECTORS - Read sens matrix by sensM/smatrix.* vectors
% S = readsensvectors(Mesh,N)

S=zeros(length(N.a),Mesh.ncells);
i=1;
for i=1:length(N.a),
    fid=fopen(['sensM' filesep 'smatrix.' num2str(i-1)],'r');
    anz=fread(fid,1,'long');
    S(i,:)=fread(fid,[1 Mesh.ncells],'double');
    fclose(fid);
end