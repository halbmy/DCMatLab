function [S,x,z]=getsens2d(datfile,x,z,N,INV)

% GETSENS2D - Loads sensitivity from matfile if possible
%             otherwise it is integrated and saved
% S = getsens2d(datfile,x,z,N);
% [S,x,z] = getsens2d(datfile);

if nargin<5, INV=[]; end
[fpath,fname,fext]=fileparts(datfile);
sensfile=strrep(datfile,fext,'-sens.mat');
issensfile=(exist(sensfile)==2);
if (~issensfile)&&(nargin<3),
    S=[];x=[];z=[];return;
end
if (nargin<2)&&(nargout>2),
    load(sensfile);
    if exist('xsave')==1, x=xsave; end
    if exist('zsave')==1, z=zsave; end
else
    calc=1;
    if issensfile,
        load(sensfile);
        %     if nargin<3,
        if exist('xsave')&&isequal(x,xsave),
            if exist('zsave')&&isequal(z,zsave),
                calc=0;
            end
        end
        %     else
        %         x=xsave;z=zsave;
        %     end
        ll=(length(x)-1)*(length(z)-1);
        if ~isequal(size(S),[length(N.r) ll]), calc=1; end
    end
    if calc,
        S=calcsens2d(x,z,N,INV);
        xsave=x;zsave=z;
        save(sensfile,'S','xsave','zsave');
    end
end