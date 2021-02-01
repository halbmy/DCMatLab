function N=read3dfile(datfile)

% READ3DFILE - Read 3d data file (all file types)
% N = read3dfile(datfile)
% the data file can be of the following types:
% 1-unified data file (see resistivity.net)
% 2-res2dinvfile
% 3-3d raw file
% 4-resecs file
% 5-pro (several 2d profiles) file (must be *.pro)
% 6-s3d (several 1d soundings) file (must be *.s3d)

N=[];
if exist(datfile)~=2, error('File does not exist!'); end
[fpath,name,ext]=fileparts(datfile);
lext=lower(ext);
if isequal(lext,'.pro'),
    N=readpro(datfile);
elseif isequal(lext,'.s3d'),
    N=readsnd3d(datfile);
else
    switch check3dfile(datfile),
        case 1,
            N=readres3dinvfile(datfile);
        case 2,
            N=readinv3dfile(datfile);
        case 3,
            N=read3drawfile(datfile);
        case 4,
            N=readresecsfile(datfile);
        case 5,
            N=readunifile(datfile);
            if isfield(N,'x')&isfield(N,'y'),
                N.elec=[N.x N.y];
                N.elec(:,3)=0;
            else
                N=readinv3dfile(datfile);
            end
        otherwise
            %error('File type unknown!');
    end 
end