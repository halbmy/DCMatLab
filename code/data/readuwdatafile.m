function DATA=readuwdatafile(datafile)

% READUWDATAFILE - Read underwater data file
% DATA = readuwdatafile(datafile)

ilon=0;ilat=0;ii=0;iu=0;ip=0;id=0;itime=0;
%%
fid=fopen(datafile);
if fid<0, 
    set(handles.statusdatafile,'String',['Could not open ' datafile]);
    return;
end
zeile=fgetl(fid);%fclose(fid);
i=0;formstr='';
while ~isempty(zeile),
    i=i+1;
    [tok,zeile]=strtok(zeile);ss='%f';
    if strcmp(lower(tok),'lon'), ilon=i;ss='%f %*s'; 
    elseif strcmp(lower(tok),'lat'), ilat=i;ss='%f %*s';
    elseif strcmp(lower(tok),'i'), ii=i;
    elseif strcmp(lower(tok),'u'), iu=i;
    elseif strcmp(lower(tok),'p'), ip=i;
    elseif strcmp(lower(tok),'d'), id=i;
    elseif strcmp(lower(tok),'channel'), ich=i;ss='%d';
    elseif strcmp(lower(tok),'time'), itime=i;ss='%s %*s';
    else ss='%*s';i=i-1; end
    formstr=[formstr ss ' '];
end
formstr(end)='';
if i==7,
    [All{1},All{2},All{3},All{4},All{5},All{6},All{7}]=textread(datafile,formstr,'headerlines',1);
elseif i==8,
    [All{1},All{2},All{3},All{4},All{5},All{6},All{7},All{8}]=textread(datafile,formstr,'headerlines',1);
else
    display('Number of fields does not match!');return;
end
% All=textscan(fid,formstr);
fclose(fid);
DATA=[];DATA.ndata=0;
if iu, DATA.u=All{iu};DATA.ndata=length(DATA.u); end
if ii, DATA.i=All{ii}; end
if id, DATA.d=All{id}; end
if ip, DATA.p=All{ip}; end
if ilat, DATA.lat=All{ilat}; end
if ilon, DATA.lon=All{ilon}; end
if itime, 
    DATA.tstr=All{itime}; 
    t0=datenum('02:00:00');
    DATA.time=datenum(DATA.tstr)+t0;
    DATA.time=DATA.time-floor(median(DATA.time));
end