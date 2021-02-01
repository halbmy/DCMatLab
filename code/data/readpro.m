function N=readpro(filename)

% READPRO - Read 3D PRO-File of 2D-Profiles
% [N,elec]=readpro(filename)
% elec..Electrode position ( x,y,z )
% N.....structure of arrays: a,b,m,n = electrode numbers(elec)
%             r = measurements    k = konfiguration factor
% PRO-File: profile1.dat x1 y1 x2 y2
%           profile2.dat x1 y1 x2 y2

input=fopen(filename,'r');
N.a=[];N.b=[];N.m=[];N.n=[];N.r=[];N.k=[];N.elec=[];
N.err=[];N.ip=[];N.i=[];N.u=[];N.rho=[];
if input<0,
    message(sprintf('Could not open profile file: %s',filename));return;
end
message(['Opening pro-file ' filename '...']);
olddir=pwd;
newdir=fileparts(filename);
sfmt='';
% if ~isempty(newdir), cd(newdir); end
try,
    profil=0;
    while 1,%count>0,
        zeile='';
        while isempty(zeile), zeile=destrip(fgetl(input));  end
        if ~ischar(zeile), break; end
        [dateiname,count]=sscanf(zeile,'%s',1);
        profil=profil+1;
        datei=[newdir filesep dateiname]; %neu
        if ~exist(datei,'file'),
            errordlg(datei,'File not found!');break;
        end
        switch check2dfile(datei),
            case 1,
                NN=readres2dinvfile(datei);
            case 2,
                NN=readinv2dfile(datei);
            case 3,
                if isempty(sfmt),
                    sfmt=inputdlg('Please specify the number of headerlines and the row numbers for A,B,M,N and the apparent resistivity (space separated)!',...
                        'Raw file format',1,{'1  1  3  5  7  12'}); 
                end
                NN=read2drawfile(datei,sfmt);
            otherwise
                errordlg(datei,'Filetype unknown!');break;
        end
        nr=length(NN.r);
        zeile(1:length(dateiname))='';
        points=sscanf(zeile,'%f %f');        
        xp=points(1:2:end);
        yp=points(2:2:end);
        if nr>0,
            message(sprintf('reading %s, (%.1f-%.1f)-(%.1f-%.1f)',datei,xp(1),yp(1),xp(end),yp(end)));
        else
            message(sprintf('2D-File %s not found!',datei));
        end
        ln=size(N.elec,1);
        % for later plotting
        %    [NN.mids,NN.seps,NN.ii,NN.kk]=midkonf2d(NN);
        % !!! problematisch wegen evtl. Löschung (ii und kk noch löschen)
        if length(xp)>length(yp), xp(end)=[]; end
        xmbm=cumsum([0;sqrt(diff(xp).^2+diff(yp).^2)]);
        elec=[interp1(xmbm,xp,NN.elec(:,1),'linear','extrap') interp1(xmbm,yp,NN.elec(:,1),'linear','extrap')];
        if size(NN.elec,2)>1, elec=[elec NN.elec(:,2)]; end
        if isfield(NN,'topo'),
            elec(:,3)=interp1(NN.topo(:,1),NN.topo(:,2),NN.elec(:,1),'linear');
            fi=find(isnan(elec(:,3)));
            if ~isempty(fi),
                last=max(find(~isnan(elec(:,3))));
                elec(fi,3)=elec(last,3);
            end
        end
        % TODO hier muss eigentlich eine abchecke hin!
        N.elec=[N.elec;elec];
        if length(unique(NN.elec(:,2)))/size(NN.elec,1)>0.3, % topo
            NN=sort2delecs(NN,1);
            di=rndig(sqrt(sum(diff(NN.elec).^2,2)),3);
            if unique(di)/size(NN.elec,1)<0.1,
               NN.elec(:,1)=[0;cumsum(di)];
               NN.elec(:,2)=0;
            end
        end
        N.zweid{profil}=NN;
        N.nr{profil}=(1:length(NN.a))+length(N.a);
        N.names{profil}=dateiname;
        N.points{profil}=points;
%         N.points{profil}=[x1 y1 x2 y2];
        %N.pro{profil}=(1:length(NN.a))+length(N.a);
        N.a=[N.a;NN.a+ln];N.m=[N.m;NN.m+ln];
        N.b=[N.b;(NN.b>0).*(NN.b+ln)];
        N.n=[N.n;(NN.n>0).*(NN.n+ln)];
        if isfield(NN,'r'), N.r=[N.r;NN.r]; end
        if isfield(NN,'rho'), N.rho=[N.rho;NN.rho]; end
        if isfield(NN,'k'), N.k=[N.k;NN.k]; end
        if isfield(NN,'err'), N.err=[N.err(:);NN.err]; end
        if isfield(NN,'ip'), N.ip=[N.ip(:);NN.ip]; end
        if isfield(NN,'i'), N.i=[N.i(:);NN.i]; end
        if isfield(NN,'u'), N.u=[N.u(:);NN.u]; end
%         [dateiname,count]=fscanf(input,'%s',1);
    end % while loop
    fclose(input);
    %     cd(olddir);
catch
    %     cd(olddir);
    fclose(input);
    display(lasterr);
    return;
end
[N.elec,I,J]=unique(round(N.elec*1000)/1000,'rows');
N.a=J(N.a);N.m=J(N.m);
fb=find(N.b);N.b(fb)=J(N.b(fb));
fn=find(N.n);N.n(fn)=J(N.n(fn));
message(sprintf('Summary of %d measurements and %d Electrodes',length(N.r),size(N.elec,1)));
if isfield(N,'err')&&(~isempty(N.err)),
    if length(N.err)==length(N.r),
        message(sprintf('Found errors in 2d files, min=%.1f%% max=%.1f%%',min(N.err)*100,max(N.err)*100));
    else
        message('Found uncomplete errors!');
    end
end
