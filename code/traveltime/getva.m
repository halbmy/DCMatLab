function [va,di]=getva(Shot)

% GETVA - Get apparent velocity of a Shot struct
% va = getva(Shot)

di=[];
if isfield(Shot,'s')&&isfield(Shot,'g'),
    if isfield(Shot,'pos'),
        di=sqrt(sum((Shot.pos(Shot.s,:)-Shot.pos(Shot.g,:)).^2,2));
    elseif isfield(Shot,'elec'),
        di=sqrt(sum((Shot.elec(Shot.s,:)-Shot.elec(Shot.g,:)).^2,2));
    else di=1;
    end
else
    for i=1:length(Shot.ns),
        rel=Shot.pos(Shot.nx{i},:)-Shot.pos(repmat(Shot.ns{i},length(Shot.nx{i}),1),:);
        di=[di;sqrt(sum(rel.^2,2))];
    end
end
va=di./Shot.t;