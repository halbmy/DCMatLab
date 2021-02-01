function [N,radius]=loadddz(infile)

% LOADDDZO - Load ddz-file
% N = loadddz(infile) or
% [N,radius] = loadddz(infile)

if nargin<1, infile='bla'; end
fid=fopen(infile,'r');
if fid<0, error('file could not be opened!'); end
for i=1:15, zeile=fgetl(fid); end
rad=str2num(zeile(findstr(zeile,':')+1:end));
N.rad=rad;
for i=1:2, zeile=fgetl(fid); end
nel=str2num(zeile(findstr(zeile,':')+1:end));
zeile=fgetl(fid);
points=str2num(zeile(findstr(zeile,':')+1:end));
for i=1:2, zeile=fgetl(fid); end
% abmniure=textscan(fid,'%d%d%d%d%f%f%f%f');
abmniure=fscanf(fid,'%f',[8 Inf])';
fclose(fid);
N.a=abmniure(:,1);
N.b=abmniure(:,2);
N.m=abmniure(:,3);
N.n=abmniure(:,4);
N.i=abmniure(:,5)/1000; %mA
N.u=abmniure(:,6)/1000; %mA
N.r=N.u./N.i; % Ohmscher Widerstand
umin=100e-6;
proz=0.03;
N.err=proz+abs(umin./N.u); % wg. mV
while length(rad)<nel, rad(end+1)=rad(1); end
for i=0:nel-1, %Süden=1, dann Uhrzeigersinn
    N.elec(i+1,1)=-rad(i+1)*sin(i/nel*2*pi);
    N.elec(i+1,2)=-rad(i+1)*cos(i/nel*2*pi);
end
N.rhoa=abmniure(:,7);
fi=find((N.err>0.15)|(N.r>-0.1));fi=[];
N.a(fi)=[];N.b(fi)=[];N.m(fi)=[];N.n(fi)=[];
N.r(fi)=[];N.u(fi)=[];N.i(fi)=[];N.err(fi)=[];
N.rhoa(fi)=[];
radius=rad(1);