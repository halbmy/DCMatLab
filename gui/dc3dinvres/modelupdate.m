function Update=modelupdate(Model,dM,modus,lbound,ubound)

% MODELUPDATE - Update Model structure Model with dM
% newModel = modelupdate(Model,dM[,modus])
% modus = 0 .. set dM
%         1 .. logarithmic update (default)
%         2 .. linear update

if nargin<5, ubound=0; end
if nargin<4, lbound=0; end
if nargin<3, modus=1; end
Update=Model;
if ~iscell(Model.M), % grid model
    switch modus,
        case 0,
            Update.M(:)=dM(:);
        case 1,
            if ubound>0,
                Update.M(:)=(ubound*(Model.M(:)-lbound).*exp(dM)+lbound*(ubound-Model.M(:)))./((Model.M(:)-lbound).*exp(dM)+ubound-Model.M(:));
            else
                Update.M(:)=(Model.M(:)-lbound).*exp(dM(:))+lbound;
            end
        case 2,
            Update.M(:)=Model.M(:)+dM(:);
    end
    return
end %else Para model
K=length(Model.z)-1;
old=0;
for k=1:K,
    sch=Update.M{k};
    kmod=prod(size(sch));
    teil=dM(old+1:old+kmod);
    switch modus,
    case 0,
        sch(:)=teil(:);
    case 1,
        if ubound>0,
            sch(:)=(ubound*(sch(:)-lbound).*exp(dM)+lbound*(ubound-sch(:)))./((sch(:)-lbound).*exp(dM)+ubound-sch(:));
        else
            sch(:)=(sch(:)-lbound).*exp(teil(:))+lbound;
        end
    case 2,
        sch(:)=sch(:)+teil(:);
    end
    Update.M{k}=sch;
    old=old+kmod;
end
