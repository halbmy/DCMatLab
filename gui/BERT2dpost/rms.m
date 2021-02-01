function erg = rms(soll,ist,lolo)

% RMS - Calculate RMS(Root Mean Square)
% rms = rms(soll,ist)
if nargin<3, lolo=0; end
% if lolo>0
%     %erg = sqrt(sum((1-log(ist)./log(soll)).^2/length(soll)))*100;
%     erg=norm(log(ist)-log(soll));
% else
%     erg = sqrt(sum((1-ist./soll).^2/length(soll)))*100;
% end
% erg=std(ist./soll)*100;
erg=sqrt(mean(((soll-ist)./soll).^2))*100;