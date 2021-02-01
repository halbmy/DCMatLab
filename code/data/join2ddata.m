function [Shot,N]=join2ddata(Shot,N)

% JOIN2DDATA - Join (2d) DC and Ra file
% [Shot,Data] = join2data(Shot,N);

pos=unique([N.elec;Shot.pos],'rows');
[tf,idc]=ismember(N.elec,pos,'rows');
[tf,itt]=ismember(Shot.pos,pos,'rows');
N.elec=pos;
N.a=idc(N.a);N.m=idc(N.m);
fi=find(N.b);N.b(fi)=idc(N.b(fi));
fi=find(N.n);N.n(fi)=idc(N.n(fi));
% showdata2d(N);
Shot.pos=pos;
if isfield(Shot,'ns'), % ns/nx style
    for i=1:length(Shot.ns), Shot.ns{i}=itt(Shot.ns{i});
        Shot.nx{i}=itt(Shot.nx{i}); end
elseif isfield(Shot,'s'), % s/g style
    Shot.g=itt(Shot.g);
    Shot.s=itt(Shot.s);
end
Shot.pos=pos;
N.elec=pos;
return
plotshot(Shot)
[pp,ff,ee]=fileparts(ttfile);
savesgtfile(fullfile(pp,['j' lower(ff) ee]),Shot);
[pp,ff,ee]=fileparts(dcfile);
saveinv2dfile(fullfile(pp,['j' lower(ff) ee]),N);

