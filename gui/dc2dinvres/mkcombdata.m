WA=readunifile('sud5m100518wa-corr.dat');
WB=readunifile('sud5m100518wb-corr.dat');
%%
WAB5=combdata2d(WA,WB);
showdata2d(WAB);%,WAB.err*100)
%%
%%
WA=readunifile('sud2m100521wa-corr.dat');
WB=readunifile('sud2m100521wb-corr.dat');
%%
WAB2=combdata2d(WA,WB);
showdata2d(WAB);%,WAB.err*100)
%%


return
%%
saveunifile('sud5m100518wabc.dat',WAB5);
saveunifile('sud2m100521wabc.dat',WAB2);