global Data DD
%%
% for i=1:length(Data),
    i=135;
    Dat=Data{i};
    showpseudoflat(Dat,abs(Dat.rho.*DD.k),20,1000);
    xlim([2 80]);
    pause(1.0);
% end