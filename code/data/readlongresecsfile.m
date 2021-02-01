function DATA = readlongresecsfile(datafile)

if nargin<1, display('Specify data file!');return; end
fid=fopen(datafile);
if fid<0, return; end
zeile=fgetl(fid);%fclose(fid);
ii=0;iu=0;ip=0;id=0;ich=1;ic1=0;ilon=0;ilat=0;itime=0;
i=0;formstr='';
while ~isempty(zeile),
    i=i+1;
    [tok,zeile]=strtok(zeile);ss='%f';
    if strcmp(lower(tok),'lon'), ilon=i;ss='%f %*s'; % N
    elseif strcmp(lower(tok),'lat'), ilat=i;ss='%f %*s';% E
    elseif strcmp(lower(tok),'i'), ii=i;
    elseif strcmp(lower(tok),'u'), iu=i;
    elseif strcmp(lower(tok),'p'), ip=i;
    elseif strcmp(lower(tok),'d'), id=i;
    elseif strcmp(lower(tok),'channel'), ich=i;ss='%d';
    elseif strcmp(lower(tok),'c1(x)'), ic1=i;ss='%f';
    elseif strcmp(lower(tok),'time'), itime=i;ss='%s %*s'; %UTC
    else ss='%*s';i=i-1; end
    formstr=[formstr ss ' '];
end
formstr(end)='';
if i>5,
    if i==6, [All{1},All{2},All{3},All{4},All{5},All{6}]=textread(datafile,formstr,'headerlines',1); end
    if i==7, [All{1},All{2},All{3},All{4},All{5},All{6},All{7}]=textread(datafile,formstr,'headerlines',1); end
    if i==8, [All{1},All{2},All{3},All{4},All{5},All{6},All{7},All{8}]=textread(datafile,formstr,'headerlines',1); end
    if i==9, [All{1},All{2},All{3},All{4},All{5},All{6},All{7},All{8},All{9}]=textread(datafile,formstr,'headerlines',1); end
else
    All=mytextscan(fid,formstr);
end
fclose(fid);
DATA=[];DATA.ndata=0;
if iu, DATA.u=All{iu};DATA.ndata=length(DATA.u); end
if ii, DATA.i=All{ii}; end
if id, DATA.d=All{id}; end
if ip, DATA.p=All{ip}; end
if ich, DATA.ch=All{ich}; end
if ic1, DATA.c1=All{ic1}; end
if ilat, DATA.lat=All{ilat}; end
if ilon, DATA.lon=All{ilon}; end
if itime, 
    DATA.tstr=All{itime}; 
    t0=datenum('00:00:00');
    DATA.time=datenum(DATA.tstr)-t0;
end 