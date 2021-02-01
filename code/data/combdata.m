function newdata=combdata(varargin)

% COMBDATA - combine two data structs (calls combdata2d/3d)
% newStruct = combdata(Struct1,Struct2);

N=varargin{1};
if isfield(N,'elec'),
  if size(N.elec,2)>2, 
    konf=combdata2d(varargin{:}); 
  else
    konf=combdata3d(varargin{:}); 
  end
else
  error('No electrodes present!');
end
