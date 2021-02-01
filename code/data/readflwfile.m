function N = readflwfile(filename)

% READ FLW (geotom flow) file
% N = readflwfile(filename)

fid=fopen(filename,'r');
if fid<0, error('file could not be opened!'); end
firstel=0;spacing=1;ippresent=0;
zeile=fgetl(fid);
while ~isempty(zeile),
    vor='';
    dp=strfind(zeile,':');
    if dp, vor=zeile(1:dp-1);nach=zeile(dp+1:end); end
    switch vor,
        case 'Name',
            name=nach;
        case 'Spacing'
            spacing=str2num(nach);
        case 'First El'
            firstel=str2num(nach);
        case 'Nr of El'
            nel=str2num(nach);
        case 'Nr of points'
            points=str2num(nach);
        case 'IP present'
            ippresent=str2num(nach);
    end
    zeile=fgetl(fid);
end
zeile=fgetl(fid);
fclose(fid);
ss='%d%d%d%d%f%f%f%f';%a b m n freq I/mA U/mV rhoa err/%
if ~isempty(strfind(zeile,';'))|~isempty(strfind(zeile,'$')), 
    ss='%d%d%d%d%*s%f%f%f%f'; end
if ippresent, 
    ss=[ss '%f']; 
    while length(sscanf(strrep(zeile,';',''),'%f'))>length(strfind(ss,'%f'))+4, ss=[ss '%f']; end
end % IP/mRad IPerr/%
fid=fopen(filename,'r');
zeile='a';while ~isempty(zeile), zeile=fgetl(fid); end
abmniure=mytextscan(fid,ss);
% abmniure=fscanf(fid,'%f',[8+ippresent*2 Inf])';
fclose(fid);
nc=length(abmniure);
ma=find(abmniure{6});
N.a=abmniure{1}(ma);%(:,1);
N.b=abmniure{2}(ma);%(:,2);
N.m=abmniure{3}(ma);%(:,3);
N.n=abmniure{4}(ma);%(:,4);
N.i=abmniure{5}(ma)/1000;%(:,5)/1000; %mA
N.u=abmniure{6}(ma)/1000;%(:,6)/1000; %mV
if nc>7, N.err=abmniure{8}(ma)/100; end %(:,6)/1000; % %
if ippresent, % convention: phases are negative
    if nc>8, N.ip=abmniure{9}(ma); end
    if nc>9, N.iperr=abmniure{10}(ma)/100; end %(:,10);  
end
iinf=65536;
N.a(N.a==iinf)=0;
N.b(N.b==iinf)=0;
N.m(N.m==iinf)=0;
N.n(N.n==iinf)=0;

nel=max([max(N.a) max(N.b) max(N.m) max(N.n)]);
N.elec=(0:nel-1)'*spacing+firstel;N.elec(:,2)=0;
N.r=abmniure{7}(ma);%(:,7);
N.k=getkonf2d(N);
% N.rho=N.u./N.i; % Ohmscher Widerstand
% umin=100e-6;
% proz=0.03;
% N.err=proz+abs(umin./N.u); % wg. mV
