function N = read3draw(datafile)
% READ2DLB - Read 3D Raw Data File
% N = read3draw('filename.dat');
% The data file must be organized as follows (column-wise)
% Ax Ay Bx By Mx My Nx Ny rho_a
% N.....structure of arrays: a,b,m,n = electrode numbers(elec)
%             r = measurements    k = konfiguration factor
%             elec = Electrode positions ( x,y,z )

%input=fopen(datafile,'r');
%if input<0, message(sprintf('Could not open datafile: %s',datafile));return; end

N.a=[];
N.b=N.a;
N.m=N.a;
N.n=N.a;
N.r=N.a;
N.elec=[];
fid=fopen(datafile,'r');
data=0;
while 1,
    zeile=fgetl(fid);
    if ~ischar(zeile), break, end
    data=data+1;
    N.b(end+1)=0;N.n(end+1)=0;
    le=size(N.elec,1);
    RR=sscanf(zeile,'%f ');
    N.r(end+1)=RR(end);RR(end)=[];
    lr2=floor(length(RR)/2);
    RRR=reshape(RR',[2 lr2])';
        [tf,loc]=ismember(RRR,N.elec,'rows');
    %tf=ismember(RRR,N.elec,'rows');
%     loc=tf*0;
%     for 1=find(loc)
%         k=0;
%         while sum(abs(N.elec(k,:)-RRR
%     end
    fi=find(~tf);
    N.elec=[N.elec;RRR(fi,:)];
    loc(fi)=(1:length(fi))+le;
    N.a(end+1)=loc(1);
    isb=(length(loc)>3);
    isn=(length(loc)>2);
    N.m(end+1)=loc(2+isb);
    if isb, N.b(end)=loc(2); end
    if isn, N.n(end)=loc(3+isb); end
end
N.elec(:,3)=0;
N.a=N.a(:);N.b=N.b(:);N.m=N.m(:);N.n=N.n(:);N.r=N.r(:);
N.k=getkonf(N);