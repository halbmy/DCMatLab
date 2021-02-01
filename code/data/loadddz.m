function [N,radius]=loadddz(infile,ccw)

% LOADDDZO - Load ddz-file
% N = loadddz(infile) or
% [N,radius] = loadddz(infile)

if nargin<2, ccw=1; end
if nargin<1, infile='bla'; end
fid=fopen(infile,'r');
if fid<0, error('file could not be opened!'); end

zeile=fgetl(fid);
while ~isempty(zeile),
    vor='';
    dp=strfind(zeile,':');
    if dp, vor=zeile(1:dp-1);nach=zeile(dp+1:end); end
    switch vor,
        case 'Radius'
            rad=str2num(nach);N.rad=rad;
        case 'First El'
            N.fel=str2num(nach);
        case 'Nr of El'
            nel=str2num(nach);
        case 'Nr of points'
            points=str2num(zeile(findstr(zeile,':')+1:end));
        case 'IP present'
            ippresent=str2num(nach);
    end
    zeile=fgetl(fid);
end
ss='%d%d%d%d%f%f%f%f';
if ippresent, ss=[ss '%f%f']; end
abmniure=mytextscan(fid,ss);
% abmniure=fscanf(fid,'%f',[8+ippresent*2 Inf])';
fclose(fid);
N.a=abmniure{1};%(:,1);
N.b=abmniure{2};%(:,2);
N.m=abmniure{3};%(:,3);
N.n=abmniure{4};%(:,4);
N.i=abmniure{5}/1000;%(:,5)/1000; %mA
N.u=abmniure{6}/1000;%(:,6)/1000; %mA
if ippresent, % convention: phases are negative
    N.ip=abmniure{8};%(:,9);
    N.iperr=abmniure{9};%(:,10); 
end
N.rho=N.u./N.i; % Ohmscher Widerstand
umin=100e-6;
proz=0.03;
N.err=proz+abs(umin./N.u); % wg. mV
while length(rad)<nel, rad(end+1)=rad(1); end
for i=0:nel-1, %Süden=1, dann Uhrzeigersinn
    N.elec(i+1,1)=-rad(i+1)*sin(i/nel*2*pi)*ccw;
    N.elec(i+1,2)=-rad(i+1)*cos(i/nel*2*pi);
end
N.r=abmniure{7};%(:,7);
% fi=find((N.err>0.15)|(N.rho>-0.1));fi=[];
% N.a(fi)=[];N.b(fi)=[];N.m(fi)=[];N.n(fi)=[];
% N.rho(fi)=[];N.u(fi)=[];N.i(fi)=[];N.err(fi)=[];
% N.r(fi)=[];
radius=rad(1);