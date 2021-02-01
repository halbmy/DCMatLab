function N=readinv3dfile(fname)

% READINV3DFILE - Read data file in inv3d format
% N = readinv3dfile(filename)
% Format: (# character can be used to comment)
% number_of_Electrodes
% x_el1 (y_el1) z_el1
% ...
% x_eln (y_eln) z_eln
% number_of_datapoints
% A_1 B_1 M_1 N_1 R_1 (Err_1)
% ...
% A_n B_n M_n N_n R_n (Err_n)
% (= Electrode numbers, 0=infinity)

N.elec=[];
fid=fopen(fname,'r');
if fid<0, error('File not found!'); end
zeile='';
while isempty(zeile), zeile=destrip(fgetl(fid)); end
ne=sscanf(zeile,'%d\n',1);
zeile=fgetl(fid);zeile(int32(zeile)==9)=' ';
if (zeile(1)=='#')&&(strfind(lower(zeile),'x')),
    ix=0;iy=0;iz=0;xmul=1;ymul=1;zmul=1;
    zeile(1)='';i=0;formstr='';
    while ~isempty(zeile),
        i=i+1;utok='';
        [tok,zeile]=strtok(zeile);
        if isempty(tok), break; end
        fis=strfind(tok,'/');
        if ~isempty(fis), % physical unit found
            fis=fis(1);utok=tok(fis+1:end);tok=tok(1:fis-1); 
        end
        fprintf('%s ',tok)
        switch lower(tok),
            case 'x', ix=i;
            case 'y', iy=i;
            case 'z', iz=i;
            case 't', it=i;
            case 'd', id=i;
        end
    end
    zeile=fgetl(fid);
end
while isempty(zeile),
  zeile=destrip(fgetl(fid));
end
for n=1:ne,
    if n>1, zeile=destrip(fgetl(fid)); end
    while isempty(zeile), zeile=destrip(fgetl(fid)); end
    el=str2num(zeile);
    N.elec(n,1:length(el))=el;
end
if size(N.elec,2)<2, N.elec(:,2)=0; end
if size(N.elec,2)<3, N.elec(:,3)=0; end
zeile='';
while isempty(zeile), zeile=destrip(fgetl(fid)); end
nm=sscanf(zeile,'%d\n',1);
zeile='';
while isempty(zeile), zeile=fgetl(fid); end
zeile(int32(zeile)==9)=' ';
if (zeile(1)=='#')&&(any(strfind(lower(zeile),'a'))),
    ia=0;ib=0;im=0;in=0;ir=0;ierr=0;iip=0;ii=0;iu=0;ik=0;it=0;
    emul=1;imul=1;umul=1;irho=0;isp=0;
    zeile(1)='';i=0;formstr='';
    while ~isempty(zeile),
        i=i+1;utok='';
        [tok,zeile]=strtok(zeile);
        fis=strfind(tok,'/');
        if ~isempty(fis), % physical unit found
            fis=fis(1);utok=tok(fis+1:end);tok=tok(1:fis-1); 
        end
        if isempty(tok), break; end
        fprintf('%s ',tok);
        switch lower(tok),
            case {'a','c1'}, ia=i;formstr=[formstr '%d'];
            case {'b','c2'}, ib=i;formstr=[formstr '%d'];
            case {'m','p1'}, im=i;formstr=[formstr '%d'];
            case {'n','p2'}, in=i;formstr=[formstr '%d'];
            case {'rhoa','ra','rho_a'}, ir=i;formstr=[formstr '%f'];
            case {'rho','r'}, irho=i;formstr=[formstr '%f'];
            case {'err','error','std'}, 
                ierr=i;formstr=[formstr '%f'];
                if isequal(utok,'%'), emul=0.01; end
            case 'ip', iip=i;formstr=[formstr '%f'];
            case 'sp', isp=i;formstr=[formstr '%f'];
            case 'i', ii=i;formstr=[formstr '%f'];
                if isequal(utok,'mA'), imul=1e-3; end
                if isequal(utok,'uA'), imul=1e-6; end
                if isequal(utok,'nA'), imul=1e-9; end
                if isequal(utok,'kA'), imul=1e+3; end
            case {'u','v'}, iu=i;formstr=[formstr '%f'];
                if isequal(utok,'mV'), umul=1e-3; end
                if isequal(utok,'uV'), umul=1e-6; end
                if isequal(utok,'nV'), umul=1e-9; end
                if isequal(utok,'kV'), umul=1e+3; end
            case 'k', ik=i;formstr=[formstr '%f'];
            case {'t','topo'}, it=i;formstr=[formstr '%f'];
            otherwise, %unknown token
                formstr=[formstr '%*s'];
        end
    end
    fprintf('found\n');
%     [ia ib im in ir ierr iip ii iu],formstr
    A=mytextscan(fid,formstr,nm,'commentstyle','#');
    if ia, N.a=A{ia}; else N.a=zeros(nm,1); end
    if ib, N.b=A{ib}; else N.b=zeros(nm,1); end
    if im, N.m=A{im}; else N.m=zeros(nm,1); end
    if in, N.n=A{in}; else N.n=zeros(nm,1); end
    if ir, N.r=A{ir}; end
    if irho, N.rho=A{irho}; end
    if iip, N.ip=A{iip}; end
    if isp, N.sp=A{isp}; end
    if ierr, N.err=A{ierr}*emul; end
    if ii, N.i=A{ii}*imul; end
    if iu, N.u=A{iu}*umul; end
    if ik, N.k=A{ik};N.konf=N.k; end
    if it, N.t=A{it}; end
else
    N.a=zeros(nm,1);N.b=N.a;N.m=N.a;N.n=N.a;N.r=N.a;N.err=N.a;N.ip=N.a;
    zeile=destrip(zeile);
    for n=1:nm,
        if n>1, zeile=destrip(fgetl(fid)); end
        while isempty(zeile), zeile=destrip(fgetl(fid)); end
        mess=str2num(zeile);
        if length(mess)<5, break; end
        N.a(n)=mess(1);N.b(n)=mess(2);N.m(n)=mess(3);N.n(n)=mess(4);
        if length(mess)>4, N.r(n)=mess(5); end
        if length(mess)>5, N.err(n)=mess(6); end
        if length(mess)>6, N.ip(n)=mess(7); end
    end
    if max(N.ip)<=0, N=rmfield(N,'ip'); end
end
fclose(fid);
nm=min([length(N.a) length(N.b) length(N.m) length(N.n)]);
fn=fieldnames(N);
for i=1:length(fn),
    fie=getfield(N,fn{i});
    if (min(size(fie))==1)&&(length(fie)>nm), fie(nm+1:end)=[];N=setfield(N,fn{i},fie); end
end
N.b(N.b<0)=0;N.n(N.n<0)=0;
% if isfield(N,'err')&&(min(N.err)<=0), 
%     messg('Found nonpositive errors. Discarding error values.');
%     N=rmfield(N,'err'); 
% end
if ~isfield(N,'k')||(length(N.k)<length(N.a)), N.k=getkonf(N); end
if ~isfield(N,'r'), % no apparent resistivities
    if ~isfield(N,'rho')&&isfield(N,'u')&&isfield(N,'i'), N.rho=N.u./N.i; end
    if isfield(N,'rho'), N.r=N.rho.*N.k; end
end
messg(sprintf('%s: %d Measurements with %d Electrodes',fname,nm,ne));

function zeile=destrip(zeile)
% strip string from comments (with # character)
aa=strfind(zeile,'#');
if ~isempty(aa), zeile=zeile(1:aa(1)-1); end
