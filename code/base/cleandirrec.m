function cleandirrec(dirname,delname)

% CLEANDIRREC - Cleans directory tree recursively
% cleandirrec(dirname,delname)
% e.g. cleandirrec('.','tmp');

if nargin<1, dirname=pwd; end
if nargin<2, delname='tmp'; end
% fprintf('%s\n',dirname);
dd=dir(dirname);
for i=1:length(dd),
    [pp,ff,ee]=fileparts(dd(i).name);
    if strcmp(ee,'.matrix')&&strcmp(dd(i).name(1),'C'),
        fn=fullfile(dirname,dd(i).name);
        fprintf('%s\n',fn); 
        delete(fn);
    end
    if (dd(i).isdir)&&~strcmp(dd(i).name,'.')&&~strcmp(dd(i).name,'..'),
        newname=[dirname filesep dd(i).name];
        if strcmp(dd(i).name,delname),
            fprintf('%s\n',newname);
            rmdir(newname,'s');
        else
            cleandirrec(newname,delname);
        end
    end
end