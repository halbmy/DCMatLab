function konf=getkonf(varargin)

% GETKONF - Get configuration factor (calls getkonf2d/3d)
% N.k = getkonf(N); % where N is data struct

N=varargin{1};
if isfield(N,'elec'),
  if size(N.elec,2)>2, 
    konf=getkonf3d(varargin{:}); 
  else
      fl=sqrt(sum(diff(N.elec([1 end],:)).^2));
      fm=sqrt(sum((N.elec(1,:)-mean(N.elec)).^2));
      if fl<fm, %apparently circle file
          konf=1;
      else
          konf=getkonf2d(varargin{:}); 
      end
  end
else
  error('No electrodes present!');
end
