function [N,ff] = readtx0file(filename,trotz)
% READTX0FILE - Read TX0 file (Lippmann multielectrode device)
% N = readtx0file(filename)
% N..structure of electrodes(elecs),numbers(a,b,m,n)
%    and possible fields (u,i,rho,r,err,ip,...)

% filename='D:\Guenther.T\2d\test060406.tx0';

if nargin<2, trotz=0; end
numel=0;del=0;firstel=0; % default values
fid=fopen(filename,'r');
if fid==-1, error('Could not open file'); end
N=[];
zeile=fgetl(fid); % read first
anfang='    A';%'#    A';
anfang2='num A';
anfang3='* num';
while ischar(zeile),
    zeile(zeile==9)=32; %tab to blank
    if ~isempty(strfind(zeile,anfang)), break; end
    if (length(zeile)>4)&&isequal(zeile(1:5),anfang2), break; end
    if (length(zeile)>4)&&isequal(zeile(1:5),anfang3), break; end
    if strfind(lower(zeile),'electrode last num'),
        numel=sscanf(zeile,'%*s%*s%*s%*s%*s%d'); 
        N.elec=(0:numel-1)'*del+firstel;N.elec(:,3)=0;
    end
    if strfind(lower(zeile),'electrode start (m)'),
        firstel=sscanf(strrep(zeile,',','.'),'%*s%*s%*s%*s%*s%f'); end
    if strfind(lower(zeile),'electrode separation (m)'),
        del=sscanf(strrep(zeile,',','.'),'%*s%*s%*s%*s%*s%f'); end
    if strfind(lower(zeile),'electrode ['),
%         aa=sscanf(strrep(zeile,',','.'),'%*s%*s%*s%*s%*s%*s%*s%*s%*s%f%f%f');
%         ii=str2num(zeile(strfind(zeile,'[')+1:strfind(zeile,']')-1));
        aa=sscanf(strrep(zeile,',','.'),'%*s%*s [%d] %*s%*s%*s%*s%*s%f%f%f');
        ii=aa(1);aa(1)=[];
        N.elec(ii,1:3)=aa;
    end
    zeile=fgetl(fid); % read next
end
if ~ischar(zeile), return; end
fr=[];
if 1,%strfind(zeile,anfang)||isequal(zeile(1:5),anfang2),
    ia=2;ib=3;im=4;in=5;ii=6;iu=7;iip=12;ierr=8;ifr=13;inr=1;
    ila=0;ilo=0;it=0;iz=0;
    zeile=destrip(zeile);
    i=0;imax=1;
    if zeile(1)=='*', zeile(1)=''; end
    while ~isempty(zeile),
        i=i+1;
        [tok,zeile]=strtok(zeile);
%         fprintf('%d %s\n',i,tok);
        switch lower(tok),
            case 'A', ia=i;
            case 'B', ib=i;
            case 'M', im=i;
            case 'N', in=i;
            case 'U', iu=i;
            case 'I', ii=i;
            case {'la','x'}, ila=i;
            case {'lo','y'}, ilo=i;
            case 'z', iz=i;
            case 'time', it=i;
            otherwise imax=i;
        end
    end
    % evaluate token string
    mai=max([ia ib im in ii iu iip ierr ila ilo iz]);
    N.a=[];N.b=[];N.m=[];N.n=[];
    if iu, N.u=[]; end
    if ierr, N.err=[]; end
    if ii, N.i=[]; end
    if ierr, N.err=[]; end
    if iip, N.ip=[]; end
    if ila, N.gpslat=[]; end
    if ilo, N.gpslon=[]; end
    if iz, N.gpsz=[]; end
    if ifr, fr=[]; end
    if inr, nr=[]; end
    if it, N.t=[]; end
    zeile=fgetl(fid);
    emul=1/100;
    imul=1/1000;umul=1/1000; %mA/mV neg.
    % evaluate unit string
    zeile=fgetl(fid);
    while ischar(zeile), % no EOF
        while ischar(zeile)&&~isempty(strfind(lower(zeile),'inf')),
            fprintf('inf ');
            zeile=fgetl(fid);
        end        
        if ~ischar(zeile), break; end
        if it==imax, 
            N.t(end+1)=datenum(zeile(max(strfind(zeile,' '))+1:end)); 
        end
        zeile(zeile==9)=32;
        bb=strrep(strrep(zeile,' - ',' 0 '),',','.');
        bb(end)=[];
        aa=str2num(bb);
%         if ~any(strfind(zeile,'- '))&&isnumeric(aa)&&(length(aa)>=mai), % valid numbers        
        if isnumeric(aa)&&(length(aa)>=mai-double(it>0))&&((aa(iu)~=0)||trotz), % valid numbers
           if ia, N.a(end+1)=aa(ia); else N.a(end+1)=0; end
           if ib, N.b(end+1)=aa(ib); else N.b(end+1)=0; end
           if im, N.m(end+1)=aa(im); else N.m(end+1)=0; end
           if in, N.n(end+1)=aa(in); else N.n(end+1)=0; end
           if ierr, N.err(end+1)=aa(ierr)*emul; else N.err(end+1)=0; end
           if iu, N.u(end+1)=aa(iu)*umul; end
           if ii, N.i(end+1)=aa(ii)*imul; end
           if iip, N.ip(end+1)=aa(iip); end % + oder - ???
           if ifr, fr(end+1)=aa(ifr); end
           if inr, nr(end+1)=aa(inr); end
           if ila, N.gpslat(end+1)=aa(ila); end
           if ilo, N.gpslon(end+1)=aa(ilo); end
           if iz&length(aa)>=iz, N.gpsz(end+1)=aa(iz); end
%            if it, N.t(end+1)=datenum(aa(inr)); end
        end
        zeile=fgetl(fid);
    end
end
fclose(fid);
% fprintf('\n');

N.ff=unique(fr);
if length(N.ff)>1, % multiple frequencies
%     N.allip=zeros(max(inr),length(N.ff));
    for i=length(N.ff):-1:1,
        fi=find(fr==N.ff(i));
        N.sip{i}=getmeasurement(N,fi);
        N.sip{i}.nr=fi;
%         N.allip(nr(fi),i)=N.ip(fi);
%         N.allrho(nr(fi),i)=N.u(fi)./N.i(fi);
    end
    aa=zeros(length(fi),1);
    N.a(nr(fi))=N.a(fi);N.b(nr(fi))=N.b(fi);
    N.m(nr(fi))=N.m(fi);N.n(nr(fi))=N.n(fi);
    N.u(nr(fi))=N.u(fi);N.i(nr(fi))=N.i(fi);
    N.ip(nr(fi))=N.ip(fi);N.err(nr(fi))=N.err(fi);
    ma=max(nr);
    N.a(ma+1:end)=[];N.b(ma+1:end)=[];N.m(ma+1:end)=[];N.n(ma+1:end)=[];
    N.u(ma+1:end)=[];N.i(ma+1:end)=[];N.ip(ma+1:end)=[];N.err(ma+1:end)=[];
end
% if iip, N.ip=atan2(N.ip,N.u); end
fn=fieldnames(N);
for i=1:length(fn),
    aa=getfield(N,fn{i});
    if min(size(aa))==1, N=setfield(N,fn{i},aa(:)); end
end
if ~isfield(N,'k')||(length(find(N.k==0))>0), N.k=getkonf2d(N); end
if ~isfield(N,'r'), % no apparent resistivities
    if ~isfield(N,'rho')&&isfield(N,'u')&&isfield(N,'i'), N.rho=N.u./N.i; end
    if isfield(N,'rho')&&isfield(N,'k')&&isequal(size(N.rho),size(N.k)),
        N.r=N.rho.*N.k; end
end
messg(sprintf('%s: %d Measurements with %d Electrodes',filename,length(N.a),size(N.elec,1)));
