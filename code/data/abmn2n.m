function N=abmn2n(AA,BB,MM,NN)

% ABMN2N - Converts position arrays for A,B,M and N into data struct
% N = abmn2n(apos,bpos,mpos,npos);
% apos etc. is array (number of data, dimension)

% [N.elec,SI,SJ]=unique([AA;BB;MM;NN],'rows');
N.elec=unique([AA;BB;MM;NN],'rows');
[TF,LOC]=ismember(AA,N.elec,'rows');
N.a=LOC;
[TF,LOC]=ismember(BB,N.elec,'rows');
N.b=LOC;
[TF,LOC]=ismember(MM,N.elec,'rows');
N.m=LOC;
[TF,LOC]=ismember(NN,N.elec,'rows');
N.n=LOC;
% [N.elec(N.a,1) N.elec(N.b,1) N.elec(N.m,1) N.elec(N.n,1)]
if size(N.elec,2)<2, N.elec(:,2)=0; end
N.k=getkonf(N);