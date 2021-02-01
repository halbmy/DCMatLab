function [S,Model]=getsens3d(datfile,Model,N,INV)

% GETSENS3D - Get 3d sensitivity from model
% either by loading it from disk or by numerical integration
% S = getsens(datfile,Model,N)
% [S,Model] = getsens(datfile)

if nargin<4, INV=[]; end
[fpath,fname,fext]=fileparts(datfile);
calc=1;
if iscell(Model.M), %Para model
    sensfile=strrep(datfile,fext,'-s.mat');
    if exist(sensfile)==2,
        load(sensfile);
        if exist('zsave')&&isequal(Model.z,zsave),
            calc=0;
        end
        calc=calc||(~isequal(size(S,1),length(N.a)));
        if calc, 
            S=modsens(Model,N); 
            zsave=Model.z;
            save(sensfile,'S','zsave');
        end
    end
else % Grid model
    sensfile=strrep(datfile,fext,'-sens.mat');
    if exist(sensfile)==2,
        xsave=[];ysave=[];zsave=[];
        load(sensfile);
        if isequal(Model.x(:),xsave(:)),
            if isequal(Model.y(:),ysave(:)),
                if isequal(Model.z(:),zsave(:)),
                    calc=0;
                end
            end
        end
        if issparse(S), calc=calc||(~isequal(size(S,2),length(N.a)));
        else calc=calc||(~isequal(size(S,1),length(N.a))); end
    end
    if calc,
        if isfield(INV,'spsens')&&(INV.spsens>0),
            if isfield(INV,'spratio'),
                S=spcalcsens3dt(Model.x,Model.y,Model.z,N,INV.spsens,INV.spratio);
            else
                S=spcalcsens3dt(Model.x,Model.y,Model.z,N,INV.spsens);
            end
        else
            S=calcsens3d(Model.x,Model.y,Model.z,N);
        end
        xsave=Model.x;ysave=Model.y;zsave=Model.z;
        save(sensfile,'S','xsave','ysave','zsave');
    end
end