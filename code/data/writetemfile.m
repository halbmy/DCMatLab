function writetemfile(filename,TEM)
% WRITETEMFILE - Write tem file for use in EM1dInv
% writetemfile(filename,TEM)

fid=fopen(filename,'w');
if isfield(TEM,'comment'), co=TEM.comment; else co='TEM file'; end
fprintf(fid,'%s\n',co);
lt=7; % loop type rectangular
fdir=3; % vertical magnetic dipole
fprintf(fid,'%d %d       !rectangular loop,vertical dipole\n',lt,fdir);
txpos=[0 0 0];rxpos=[0 0 0]; % Tx pos
if isfield(TEM,'txpos'), txpos=TEM.txpos; end
if isfield(TEM,'rxpos'), rxpos=TEM.rxpos; end
fprintf(fid,'%g %g %g %g %g %g    !Tx position Rx position\n',txpos,rxpos);
fprintf(fid,'%d %d   ! loop dimension\n',TEM.loopsize);
fprintf(fid,'%d %d %d   ! dB/dt for input, output and inversion\n',3,3,3);
wft=0;ntwf=0; % step
fprintf(fid,'%d %d  ! waveform (0-step,1-impulse,3-user-def.\n',wft,ntwf);
if wft>1, % step
    for i=1:ntwf, fprintf(fid,'...'); end    
end
fprintf(fid,'%d %d 0 ! no filters, no frontgate\n',0,0);
A=zeros(TEM.ndata,5);
A(:,1)=TEM.t;
% A(:,2)=TEM.v;
A(:,2)=TEM.v/TEM.coilsize; % dB/dt
if isfield(TEM,'err'), A(:,3)=TEM.err; else A(:,3)=0.05; end
fprintf(fid,'%g\t%g\t%g\t%d\t%d\n',A');
fclose(fid);