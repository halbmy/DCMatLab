function N = readsyscalfile(filename)

% READ IRIS Syscal Pro file
% Data = readsyscalfile(filename)
% Author: Tobias Pfaff (Uni Heidelberg)

fin = fopen(filename,'r');

% check header
id = fread(fin,20,'uchar')';
cid = [2 0 0 128 8 'SYSCAL Pro' 0 ' at '];
if (id~=cid); error ('not a valid syscal pro file'); end;

% enter main part
fseek(fin,1029,'bof');

elec=[]; abmn=[]; data=[];
while (1)
    % read block        
    [atype,cnt] = fread(fin,1,'uint16'); % array type : 0 W 4 WS 6 DD 
    if (cnt == 0) break; end;
    aa=fread(fin,1,'uint16');
%     if(fread(fin,1,'uint16') ~= 3);error('unknown parameter');end;
    if ~ismember(aa,[2 3]), display('unknown parameter');break; end
    mtime = fread(fin,1,'float32'); % measurement time [ms]
    m_delay = fread(fin,1,'float32'); % delay before chargeability measurement
    bb=fread(fin,1,'uint32');
    if( bb ~= 1);error('unknown parameter');end;

    % find electrode positions in table
    x_pos = fread(fin,4,'float32')'; % c1 c2 p1 p2
    y_pos = fread(fin,4,'float32')'; % c1 c2 p1 p2
    z_pos = fread(fin,4,'float32')'; % c1 c2 p1 p2
    if (any(y_pos))
        break;
    end
    for n=1:4
        if (isempty(elec) || ~any(elec(:,1)==x_pos(n) & (elec(:,2)==y_pos(n))))
            elec=[elec;x_pos(n) y_pos(n)];
            e_num(n) = size(elec,1);
        else
            e_num(n) = find(elec(:,1)==x_pos(n) & (elec(:,2)==y_pos(n)));
        end
    end
    abmn = [abmn; e_num];

    sp = fread(fin,1,'float32'); % spontanious polarisation
    vp = fread(fin,1,'float32'); % voltage difference
    in = fread(fin,1,'float32'); % injected current
    rho= fread(fin,1,'float32'); % resisitivity
    gm = fread(fin,1,'float32'); % global chargeability
    dev= fread(fin,1,'float32') * 0.01; % std. deviation    
    t_m = fread(fin,20,'float32'); % chareability time windows
    m_x = fread(fin,20,'float32'); % partial chargeability    
    status = fread(fin,1,'uint16'); % status bits (80:multichannel(lower bits=channel #) 16:single channel)
    num = fread(fin,1,'uint16'); % number of measurement (starting with 0)
    name = fread(fin,20,'schar')'; name(name==0)=[];    
    stacks = fread(fin,1,'float32'); % number of stacks measured
    rs_check = fread(fin,1,'float32'); % rs_check reception dipole
    vab = fread(fin,1,'float32'); % absolute injected voltage
    bat_tx = fread(fin,1,'float32'); % tx battery voltage
    bat_rx = fread(fin,1,'float32'); % rx battery voltage
    temp = fread(fin,1,'float32'); % temperature
    if aa==3, dtime = fread(fin,2,'uint32'); end % date & time in some strange format    
    %gtime = (dtime(2)*(2^32/fact)+dtime(1)/fact - d0)/24/3600 - datenum(2004,0,0); % date in day-since-04 format
    data = [data; rho dev in vp gm];
end;
fclose(fin);

% resort electrodes
[selec,ie] = sortrows(elec);
sabmn = abmn;
for n=1:length(ie)
    sabmn(abmn(:)==ie(n)) = n;
end

N = struct('a',sabmn(:,1),'b',sabmn(:,2),'m',sabmn(:,3),'n',sabmn(:,4), ...
           'elec',selec,'r', data(:,1), 'err', data(:,2), 'i', data(:,3)/1000, ...
           'u', data(:,4)/1000, 'ip', data(:,5));
N.k=getkonf2d(N);
N.r=N.u ./ N.i .* N.k;