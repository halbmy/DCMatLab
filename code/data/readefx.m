function N = readefx(datfile)

% READEFX - Read geosys device EFX file
% N = readefx(filename)

fid=fopen(datfile,'r');
zeile=fgetl(fid);
while isempty(zeile)||(zeile(1)==39), %' character
    zeile=fgetl(fid);
end
N.a=[];N.b=[];N.m=[];N.n=[];N.i=[];N.u=[];N.r=[];N.ip=[];
ia=1;im=2;in=3;ib=4;ii=10;iu=11;ir=12;iip=13;
zeile=fgetl(fid); % token row !!! has to be interpreted
while ischar(zeile),
   z1=strrep(strrep(zeile,';;',';0;'),';;',';');
   aa=str2num(strrep(strrep(strrep(z1,',','.'),':',''),';',' '));
   if ia, N.a(end+1)=aa(ia); end
   if ib, N.b(end+1)=aa(ib); end
   if im, N.m(end+1)=aa(im); end
   if in, N.n(end+1)=aa(in); end
   if ii, N.i(end+1)=aa(ii)/1000; end
   if iu, N.u(end+1)=aa(iu)/1000; end
   if ir, N.r(end+1)=aa(ir); end
   if iip, N.ip(end+1)=aa(iip); end
   zeile=fgetl(fid);
end
N.a=N.a(:);N.b=N.b(:);N.m=N.m(:);N.n=N.n(:);
N.i=N.i(:);N.u=N.u(:);N.r=N.r(:);N.ip=N.ip(:);
mm=max([N.a;N.b;N.m;N.n]);
fclose(fid);
N.elec=(0:mm-1)'*0.05;N.elec(:,2)=0;
N.k=getkonf2d(N);
% rhoa=N.u./N.i.*N.k;
% plot(N.r,rhoa,'.')
% showdata2d(N,N.r,struct('clog',1))