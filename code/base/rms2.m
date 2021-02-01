function erg = rms(soll,ist,lolo)

% RMS - Calculate RMS(Root Mean Square)
% rms = rms(soll,ist,lolo)
% lolo - logarithmic
erg=std(ist./soll-1)*100;
% if nargin<3, lolo=0; end
% if lolo>1,
%     %erg = sqrt(sum((1-log(ist)./log(soll)).^2/length(soll)))*100;
%     erg=norm(log(ist)-log(soll))/length(soll)*100;
% else
%     erg = sqrt(sum((1-ist./soll).^2/length(soll)))*100;
% end