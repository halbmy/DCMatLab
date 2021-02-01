function [R,Rez,MEA]=mfdfwd2d(Mod,N,FOR)

% FDFWD2D - 2D DC Forward Calculation with finite differences
% Rhoa = mfdfwd3d(Model,N,OPT)
% Model    - model structure containing x,z (grid nodes) and M (model)
%   N      - Structure of electrode numbers(a,b,m,n), 
%            k-factors(k) and measurements(r)
%            elec- Electrode Positions
%   OPT    - structure of possible fields
%            refine - 
%            acc - Accuracy of forward step
%            tol - Tolerance for incomplete preconditioner
%            maxit - Maximum Iterations for pcg
%            rand - boundaring cells (3)
%            prolong - prolonging factor (5)
%            zusatz - cells outside Electrodes (2)

if nargin<2, error('At least 2 input arguments!'); end
if nargin<3,
    FOR=struct('method',0,'acc',1e-3,'tol',1e-4,'maxit',50,...
        'rand',4,'prolong',5,'zusatz',4,'direct',-1); 
end
Lay=0;
if isfield(Mod,'Lay'), Lay=Mod.Lay; end
if nargout>2, [R,Rrez,MEA]=dcfwd2d(Mod.x,Mod.z,Mod.M,Lay,N,FOR); 
elseif nargout>1, [R,Rrez]=dcfwd2d(Mod.x,Mod.z,Mod.M,Lay,N,FOR);
else R=dcfwd2d(Mod.x,Mod.z,Mod.M,Lay,N,FOR); end