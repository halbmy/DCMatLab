function [tauopt,appnew]=linesearch(soll,old,new)

% LINESEARCH - Applies line-search procedure
% tauopt = linesearch(old,new)

for i=1:40,
    tau=i*0.05;
    inbet=old+tau*(new-old);
    ch(i)=sqrt(mean((soll-inbet).^2));
%     ch(i)=rms(soll,inbet);
end
[xx,nn]=min(ch);
tauopt=nn*0.05;    
appnew=old+tauopt*(new-old);