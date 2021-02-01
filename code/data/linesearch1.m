function [tauopt,appR]=linesearch(N,oldR,R,islog)

% LINESEARCH - Applies line-search procedure
% tauopt = linesearch(N,oldR,newR)
% [tauopt,appR] = ...

if nargin<4, islog=0; end
for i=1:20,
    tau=i*0.05;
    if islog, appR=oldR.*exp(tau*(log(R)-log(oldR))); 
    else appR=oldR+tau*(R-oldR); end
    ch(i)=chi2(N.r,appR,N.err,1);
    %if nargin>3, %+ model functional
end
[xx,nn]=min(ch);
tauopt=nn*0.05;    
appR=oldR+tauopt*(R-oldR);