function [S,x,y,z]=get3dsens(datfile,N,x,y,z)

% GET3DSENS - Loads sensitivity from matfile if possible
%           otherwise it is integrated and saved
% S = get3dsens(datfile,N,x,y,z);
% [S,x,y,z] = get3dsens(datfile);

[fpath,fname,fext]=fileparts(datfile);
sensfile=strrep(datfile,fext,'-sens.mat');
calc=1;
if (exist(sensfile)==2),
    load(sensfile);
    if nargin>3,
        if exist('xsave')&&isequal(x,xsave),
            if exist('ysave')&&isequal(y,ysave),
                if exist('zsave')&&isequal(z,zsave),
                    calc=0;
                end
            end
        end
    else
        x=xsave;y=ysave;z=zsave;
        ll=(length(x)-1)*(length(z)-1);
        calc=~isequal(sort(size(S)),sort([length(N.r) ll]));
    end
end
if (nargin>3)&calc,
    S=calcsens3d(x,y,z,N);
    xsave=x;ysave=y;zsave=z;
    save(sensfile,'S','xsave','ysave','zsave');
end
