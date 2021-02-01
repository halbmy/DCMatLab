function [ELE,NODE,FACE]=loadmesh2d(meshname)

% LOADMESH - Load DCFEMLIB 2d ASCII mesh
% [ELE,NODE,FACE]=loadmesh3d(meshname)

if 1,
if nargin<1, error('Meshname must be specified'); end
ca=exist([meshname '.e']); % carsten-format (ohne 1. Zeile und Spalte)
if ca,
    fid=fopen([meshname '.e'],'r');
    zeile=fgetl(fid);
    fclose(fid);
    spalten=length(str2num(zeile));
    fid=fopen([meshname '.e'],'r');
    ELE=fscanf(fid,'%d',[spalten Inf])'+1;
    fclose(fid);
    if spalten==8, ELE=ELE(:,2:4); end
%         ELE(:,1)=(1:size(ELE,1))'; end
    fid=fopen([meshname '.n'],'r');
    NODE=fscanf(fid,'%f',[3 Inf])';
    %NODE(:,3)=abs(NODE(:,3));
    fclose(fid);
else
    fid=fopen([meshname '.ele'],'r');
    si=fscanf(fid,'%d',[3 1]);
    ELE=fscanf(fid,'%d',[6 si(1)])'+1;ELE(:,1)=[];
    fclose(fid);
    fid=fopen([meshname '.node'],'r');
    si=fscanf(fid,'%d',[4 1]);
    NODE=fscanf(fid,'%f',[5 si(1)])';
    NODE(:,1)=[];%NODE(:,3)=abs(NODE(:,3));
    fclose(fid);
end
if nargout>2,
    if ca,
        fid=fopen([meshname '.f'],'r');
        FACE=fscanf(fid,'%f',[6 Inf])';
        fclose(fid);
        FACE(:,4:5)=[];
    else
        fid=fopen([meshname '.face'],'r');
        si=fscanf(fid,'%d',[2 1]);
        FACE=fscanf(fid,'%f',[5 si(1)])';
        fclose(fid);
        FACE(:,1)=[];
    end
    fi=find(FACE(:,4)==-1); 
    %==-1 for neumann, ==-2 for dirichlet, <0 for both
    FACE=FACE(fi,1:3)+1;
end
else
if nargin<1, error('Meshname must be specified'); end
ca=exist([meshname '.e']); % carsten-format (ohne 1. Zeile und Spalte)
if ca,
    fid=fopen([meshname '.e'],'r');
    ELE=fscanf(fid,'%d',[5 Inf])'+1;
    fclose(fid);
    fid=fopen([meshname '.n'],'r');
    NODE=fscanf(fid,'%f',[4 Inf])';
    NODE(:,3)=abs(NODE(:,3));
    fclose(fid);
else
    fid=fopen([meshname '.ele'],'r');
    si=fscanf(fid,'%d',[3 1]);
    ELE=fscanf(fid,'%d',[6 si(1)])'+1;ELE(:,1)=[];
    fclose(fid);
    fid=fopen([meshname '.node'],'r');
    si=fscanf(fid,'%d',[4 1]);
    NODE=fscanf(fid,'%f',[5 si(1)])';
    NODE(:,1)=[];NODE(:,3)=abs(NODE(:,3));
    fclose(fid);
end
if nargout>2,
    if ca,
        fid=fopen([meshname '.f'],'r');
        FACE=fscanf(fid,'%f',[6 Inf])';
        fclose(fid);
        FACE(:,4:5)=[];
    else
        fid=fopen([meshname '.face'],'r');
        si=fscanf(fid,'%d',[2 1]);
        FACE=fscanf(fid,'%f',[5 si(1)])';
        fclose(fid);
        FACE(:,1)=[];
    end
    fi=find(FACE(:,4)==-1); 
    %==-1 for neumann, ==-2 for dirichlet, <0 for both
    FACE=FACE(fi,1:3)+1;
end
end