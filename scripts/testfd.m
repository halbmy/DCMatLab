x=-7:7;
y=-5:5;
z=0:6;
SIGMA=ones(length(x)+1,length(y)+1,length(z)+1);
SIGMA(:,:,1)=0;
C1=diskr_dm(x,y,z,SIGMA);
taucs_ccs_write_binary(C1,'C1.bin');
SIGMA=SIGMA/100;
SIGMA(3:5,3:5,3:5)=SIGMA(3:5,3:5,3:5)*2;
C=diskr_dm(x,y,z,SIGMA);
taucs_ccs_write_binary(C,'C.bin');
ELEC=[0 0 0];
MEA=fd3dmea(x,y,z,SIGMA,ELEC,C);
