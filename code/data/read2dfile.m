function N=read2dfile(datfile)

% READ2DFILE - Read 2d data file (all file types)
% N = read2dfile(datfile)

N=[];
if exist(datfile)~=2, error('File does not exist!'); end
[fpath,name,ext]=fileparts(datfile);
lext=lower(ext);
if isequal(lext,'.amp'), % ABEM multi-purpose file
    N=readampfile(datfile);
elseif isequal(lext,'.stg'), % STING device file
    N=readstgfile(datfile);
elseif isequal(lext,'.flw'), % GEOTOM device file
    N=readflwfile(datfile);
elseif isequal(lext,'.tx0'), % Lippmann device file
    N=readtx0file(datfile);
elseif isequal(lext,'.bin'), % IRIS SYSCAL Pro device file
    N=readsyscalfile(datfile);
elseif isequal(lext,'.s2d'), % SOUNDINGS file
    N=readsnd2d(datfile);
elseif ismember(lext,{'.wen','.dd','.pd','.pp','.slm'})
    N=readres2dinvfile(datfile);
else
    cf=check2dfile(datfile); %1-res2dinv,2-unified data,3-raw,4-resecs
    % if unknown & *.txt -> resecs
    if (cf==0)&&isequal(lower(ext),'.txt'), cf=4; end %resecs
%     if findex==4, cf=4; end %does not work on v6.5
    switch cf,
        case 1, % RES2dINV file
            N=readres2dinvfile(datfile);
        case 2, % house format
            N=readinv2dfile(datfile);
%             N=readunifile(datfile);
            if isfield(N,'i')&&(max(N.i)>10), N.i=N.i/1000; end %A statt mA
            if isfield(N,'u')&&(max(N.u)>10), N.u=N.u/1000; end %V statt mV
        case 3, % raw data file
            sfmt=inputdlg('Please specify the number of headerlines and the row numbers for A,B,M,N and the apparent resistivity (space separated)!',...
                'Raw file format',1,{'1  1  4  7  10  13'});
            if ~isempty(sfmt), N=read2drawfile(datfile,str2double(sfmt{1})); end
        case 4, % RESECS device file
            N=readresecsfile(datfile);
            if isfield(N,'u'), fi=find(abs(N.u)>9.9); else fi=[]; end
            if ~isempty(fi),
                fn=fieldnames(N);
                for i=1:length(fn),
                    fie=getfield(N,fn{i});
                    if size(fie,2)==1, fie(fi)=[]; end
                    N=setfield(N,fn{i},fie);
                end
            end
            NN=getpseudos(N);
            nr=1;le=length(NN.zweid);
            if le==0, display('No profiles found!');return; end
            if le>1,
                ss=inputdlg([num2str(le) 'profiles present. Specify number!']);
                nr=round(str2double(ss{1}));
            end
            N=NN.zweid{nr};
        otherwise
            display('File format unknown!');
            return;
    end
end
fi=find(N.a==0);
if any(fi), 
    N.a(fi)=N.b(fi);N.b(fi)=0;
    if isfield(N,'k'), N.k(fi)=-N.k(fi); end
    if isfield(N,'k'), N.u(fi)=-N.u(fi); end
end

% switch check2dfile(datfile),
%     case 1,
%         N=readres2dinvfile(datfile);
%     case 2,
%         N=readinv2dfile(datfile);
%     case 3,
%         sfmt=inputdlg('Please specify the number of headerlines and the row numbers for A,B,M,N and the apparent resistivity (space separated)!',...
%             'Raw file format',1,{'1  1  3  5  7  12'});
%         N=read2drawfile(datfile,sfmt);
%     case 4,
%         NN=readresecsfile(datfile);
%         NN=getpseudos(NN);nr=1;le=length(NN.zweid);
%         if le==0, display('No profiles found!');return; end
%         if le>1, 
%             ss=inputdlg([num2str(le) 'profiles present. Specify number!']);
%             nr=str2double(ss{1});
%         end
%         N=NN.zweid{nr};
%     otherwise
%         error('File type unknown!');
% end 