function [err,out]=estimateerror(N,proz,Umin,I)

% ESTIMATEERROR - Estimate Errors (and noisyfies data)
% error = estimateerror(N,proz,Umin,I)
% [N.err,noisy_data] = estimateerror(N,proz,Umin,I)
% N - data structure with data(N.r) and konf. factors(N.k)
% proz - relative error in percent (1%) - scalar or vector 
% Umin - minium voltage (100µV)
% I - driving current (100mA)

if nargin<1, error('Data set required'); end
if nargin<2, proz=3; end
if nargin<3, Umin=100e-6; end
if nargin<4, I=100e-3; end

if isfield(N,'i')&&(length(N.i)==length(N.k)), I=N.i; end
if isfield(N,'u')&&(length(N.u)==length(N.k)),
    err=abs(Umin./N.u)+proz/100;
else
    if isfield(N,'rho')&&(length(N.rho)==length(N.a)),
        err=abs(1./N.rho)*Umin./I+proz/100;
    elseif isfield(N,'r')&&(length(N.r)==length(N.a)),
        err=abs(N.k./N.r)*Umin./I+proz/100;
    else
        display('Error estimation: Could not get voltage information!');
        err=ones(size(N.a))*proz/100;
    end
end
if nargout>1,
  noise=randn(size(N.a)).*err;  
  if isfield(N,'r'), out=N.r.*(1+noise); end
  message(sprintf('Error min=%.1f%% max=%.1f%% mean=%.1f%%',...
      min(abs(noise))*100,max(abs(noise))*100,mean(abs(noise))*100));
end
