function [tauopt,appR]=linesearch(N,oldR,R,islog,num)

% LINESEARCH - Applies line-search procedure
% tauopt = linesearch(N,oldR,newR)
% [tauopt,appR] = ...

if nargin<5, num=20; end
if nargin<4, islog=1; end
ch0=chi2(N.r,oldR,N.err,islog);
di=1/num;
ch=zeros(num,1);
for i=1:num,
    tau=i*di;
    if islog, appR=oldR.*exp(tau*(log(R)-log(oldR))); else appR=oldR+tau*(R-oldR); end
    ch(i)=chi2(N.r,appR,N.err,islog);
    %if nargin>3, %+ model functional
end
[xx,nn]=min(ch);
tauopt=nn*di;
if ch(nn)>ch0, tauopt=0; end
if islog, appR=oldR.*exp(tauopt*(log(R)-log(oldR))); else appR=oldR+tauopt*(R-oldR); end
