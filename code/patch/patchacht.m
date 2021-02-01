function patchacht(x,y,z,fcol)

if nargin<4, fcol='b'; end
Faces = [ ...
  1 2 6 5 ;
  2 3 7 6 ;
  3 4 8 7 ;
  4 1 5 8 ;
  1 2 3 4 ;
  5 6 7 8 ];
Vert=[x(:) y(:) z(:)];
ptch=patch('Vertices',Vert,'Faces',Faces,'FaceColor','b'); 
