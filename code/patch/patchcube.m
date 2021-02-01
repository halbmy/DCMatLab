function ptch = patchcube(x1,x2,y1,y2,z1,z2,fcol)

% PATCHCUBE - Patch cube
% ptch = patchcube(x1,x2,y1,y2,z1,z2[,facecolor]

if nargin<7, fcol=[1 0.5 0.5]; end

Vert = [ ...
  x1 y1 z1 ;
  x2 y1 z1 ;
  x2 y2 z1 ;
  x1 y2 z1 ;
  x1 y1 z2 ;
  x2 y1 z2 ;
  x2 y2 z2 ;
  x1 y2 z2]; 
Faces = [ ...
  1 2 6 5 ;
  2 3 7 6 ;
  3 4 8 7 ;
  4 1 5 8 ;
  1 2 3 4 ;
  5 6 7 8 ];
ptch=patch('Vertices',Vert,'Faces',Faces,...
    'FaceVertexCData',jet(64),'FaceColor',fcol); 
