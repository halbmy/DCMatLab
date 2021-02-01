function wc=veltocrim(velocity,ematrix,n)

% VELTOCRIM - Compute water content from velocity by CRIM formula
%             (valid for porous media without clay)
% veltocrim(velocity,ematrix,porosity)
% water_content = (v_0/v + (n-1)*sqrt(eps_matrix)-n)/(sqrt(eps_water)-1)
% where v..velocity, n..porosity, eps_matrix..matrix permittivity
% wc = veltocrim(velocity,eps_matrix,porosity)
%      the latter two can be neglected (default: eps_matrix=5,n=40%)

if nargin<1, error('Specify velocity!'); end
if nargin<3, n=0.4; end
if nargin<2, ematrix=5; end
v0=3e8;ewater=81;
wc=(v0./velocity+(n-1)*sqrt(ematrix)-n)/(sqrt(ewater)-1);
