function [Alldata,freq]=readradicresfile(resfile,rememcoup,fcoupl)

% READRADICRESFILE - Read *.RES file from Radic Research SIP256C device
% [Alldata,freq] = readradicresfile(filename[,removeemcoupling])

if nargin<3, fcoupl=1:3; end
if nargin<2, rememcoup=0; end

Data=struct('a',[],'b',[],'m',[],'n',[],'elec',[]);
fid=fopen(resfile); % open file for first time
if fid==-1, error('Could not open file!');
end
zeile='';
% read number of remote units
while isempty(findstr(zeile,'Number of Remote Units')),
    zeile=fgetl(fid);
end
nru=str2num(zeile(find(zeile==']')+1:end));
while isempty(findstr(zeile,'Begin Layout')),
    zeile=fgetl(fid);
end
% reading the electrode positions from the file
for i=1:nru,
    aa=str2num(fgetl(fid));
    Data.elec(i,:)=aa(3:4);
end
fclose(fid);
%% read array configuration
Data.a=[];Data.b=[];Data.m=[];Data.n=[];
fid=fopen(resfile); % reopen file (readings come first)
zeile='';
while isempty(findstr(zeile,'Number of Readings')),
    zeile=fgetl(fid);
end
% get number of readings from file
nread=str2num(zeile(find(zeile==']')+1:end));
zeile=fgetl(fid); % begin readings
lnum=0;
for i=1:nread,
    zeile=fgetl(fid); %read line
    pos=str2num(zeile);pos(1)=[]; % numbers, first is equal to i
    gegen=(1:nru)'+1; % electrode to measure against for each el.
    for n=3:length(pos), % all potential bridges
        fi=find(gegen==pos(n)); % against electrode = bridge
        gegen(fi)=gegen(fi)+1; % increase by one
    end
    gegen(pos(1:2))=0; % current electrodes cannot be potential
    gegen(max(pos(1:2)-1,1))=0; % current-1 cannot be potential
    POS{i}=double(gegen>0);
    lfi=sum(POS{i});
    POS{i}(find(gegen>0))=(1:lfi)'+lnum;
    lnum=lnum+lfi;
    fi=find(gegen); % find only valid against electrodes
    Data.a=[Data.a;ones(length(fi),1)*pos(1)]; % constant A
    Data.b=[Data.b;ones(length(fi),1)*pos(2)]; % constant B
    Data.m=[Data.m;fi];        % M is position (index) itself
    Data.n=[Data.n;gegen(fi)]; % N is against electrode
end
while isempty(findstr(zeile,'Begin FrequencyParameter')), zeile=fgetl(fid); end
freq=[];
zeile=fgetl(fid);
while isempty(findstr(zeile,'End FrequencyParameter')),
    fp=str2num(zeile);
    if (length(fp)<7)||(fp(7)>0), freq(end+1)=fp(1); end
    zeile=fgetl(fid);
end
fclose(fid);
%% read actual data
fid=fopen(resfile); % reopen file (readings come first)
Data.r=ones(size(Data.a));Data.ip=Data.r;Data.err=Data.r;
Data.iperr=Data.r;
Alldata={};for i=1:length(freq), Alldata{i}=Data; end
reading=0;runit=0;pp=0;
ipvals=zeros(size(freq));
iperr=zeros(size(ipvals));
opt=struct('fdamp',1e3);
while 1,
    zeile=fgetl(fid);
    while ~isempty(zeile)&&(zeile(1)==' '), zeile(1)=''; end
    if ~ischar(zeile), break; end
    if (length(zeile)>7)&&isequal(zeile(1:8),'Reading:'),
        reading=str2num(zeile(9:14)); end
    if (length(zeile)>11)&&isequal(zeile(1:12),'Remote Unit:'),
        runit=str2num(zeile(13:end)); end
    if (length(zeile)>4)&&isequal(zeile(1:5),'Freq.'),%
        if (reading>0)&(runit>0), pp=POS{reading}(runit); end
        if (reading>0)&&(runit>0)&&(pp>0),
            for i=1:length(freq),
                zeile=fgetl(fid);ii=min(findstr(zeile,':'));
                if ~isempty(ii)&&(ii>1), zeile=zeile(1:ii-1); end
                if (length(zeile)>35)&&(zeile(36)~=' '), zeile(35)=' '; end
                if (length(zeile)>25)&&(zeile(26)~=' '), zeile(25)=' '; end
                nn=str2num(strrep(strrep(zeile,'c',''),'n','')); % ,' . ',''
                if isempty(nn), break; end
                ipvals(i)=nn(3)*1000/180*pi;
                iperr(i)=nn(5)*1000/180*pi;
                if length(nn)>4,
                    Alldata{i}.r(pp)=nn(2);
                    Alldata{i}.ip(pp)=-nn(3)*1000/180*pi; %mrad positive
                    Alldata{i}.err(pp)=nn(4)/100;
                    Alldata{i}.iperr(pp)=nn(5)*1000/180*pi; %mrad
                end %enough values
                if length(nn)>6,
                    Alldata{i}.konf(pp)=nn(7);
                    Alldata{i}.imp(pp)=nn(2)/nn(7);
                else
                    a=1;
                end
            end %for i=1:length(freq)
            % do the EM coupling
            %%
            if rememcoup>0,
                ipcorr=ipvals;ii=find(ipvals~=0);
                ipcorr(ii)=removeemcoupling(freq(ii),-ipvals(ii)/1000,iperr(ii)/1000,opt)*1000;
                %                 semilogx(freq,-ipvals,freq,ipcorr);
                for i=1:length(freq), Alldata{i}.ip(pp)=ipcorr(i); end
                %                removeemfromphase(freq,ipvals,fcoupl,removeemcoupling,1e99);pause(0.1);
                %                [ipcorr,ff,ascent]=removeemfromphase(freq,-ipvals,fcoupl,removeemcoupling,1e99);
                %                if ascent>0.5, % successfull
                %                    for i=1:length(freq), Alldata{i}.ip(pp)=-ipcorr(i); end
                %                    nemremoved=nemremoved+1;
                %                end
            end
        end
    end
end
fclose(fid);