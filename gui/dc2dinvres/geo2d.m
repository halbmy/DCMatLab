% phi = geo2d(x,z,sigma,sx,sz,px,pz)
%
% User interface for mex-routine geo2dmx by C. Schwarzbach and R.-U. Boerner.
% Calculates the electrical potential phi for a number of source Electrodes ns
% given by vectors sx and sz and potential Electrodes given by vectors px and pz. 
% The model is provided by vectors x, z of lengths nx, nz which contain the 
% coordinates for a rectangular grid. (nz-1)x(nx-1) matrix sigma contains
% the cell conductivities for this grid.
% (np)x(ns) matrix phi returns potentials at all Electrodes, the is-th column
% of phi corresponding to the is-th source position, the ip-th row of phi
% corresponding to the ip-th receiver location.

function phi = geo2d(x,z,sigma,sx,sz,px,pz)
% set default values for control paremeters
if nargin<4, error('At least 4 input arguments required!'); end
if nargin<5, sz=sx*0; end
if nargin<6, px=sx;pz=sz; end
if nargin==6, pz=px*0; end
NoProlong = 0;    % No. of grid nodes to be added at the boundaries
FacProlong = 5.0; % Factor to stretch additional cells at the boundary
ModProlong = 'b'; % Switch for choice of prolonging conductivities ([a]verage or [b]oundary)
NoRefine = 0;     % No. of grid nodes to refine each cell row and column
NoLegendre = 8;   % No. of Gauss-Legendre abscissas
NoLaguerre = 8;   % No. of Gauss-Laguerre abscissas
Boundary = 'm';   % Switch for boundary condition ([d]irichlet or [m]ixed)
Background = 's'; % Choice of background conductivity for singularity removal:
                  % 'a' ... average conductivity of the whole (nonprolonged) grid
                  % 'm' ... average (mean) conductivity at the source location
                  %         (maximum of four different conductivities averaged 
                  %          if source is located on a node)
                  % 's' ... conductivity of source location - exact singularity removal

phi = geo2dmx(x,z,sigma,sx,sz,px,pz);
%phi = geo2dmx(x,z,sigma,sx,sz,px,pz, ...
%               NoProlong,FacProlong,ModProlong,NoRefine, ...
%               NoLegendre,NoLaguerre,Boundary,Background);
