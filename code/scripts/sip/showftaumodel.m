function showftaumodel(A,z,tau)

% SHOWFTAUMODEL - Show SIP f-tau model as imageplot
% showftaumodel(A,z,tau)

imagesc(A);colorbar;
set(gca,'XTick',1:2:length(tau),'XTickLabel',num2strcell(rndig(tau(1:2:end)*1000,2)));
set(gca,'YTick',0.5:length(z)-0.5,'YTickLabel',num2strcell(z));
xlabel('\tau in ms');
ylabel('depth in m');
title('spectral chargeability');
