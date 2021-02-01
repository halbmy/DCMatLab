if 0,
    dcfile='d:\Guenther.T\2d\parkdeck\Parkdeck2.dat';
    ttfile='d:\Guenther.T\2d\parkdeck\parkdeck.sgt';
else
    [fname,pname]=uigetfile('*.dat;*.ohm');if ~ischar(fname), return; end
    dcfile=fullfile(pname,fname);
    [fname,pname]=uigetfile('*.dat;*.sgt');if ~ischar(fname), return; end
    ttfile=fullfile(pname,fname);
end
N=readinv2dfile(dcfile);
Shot=readunishot(ttfile);
pos=unique([N.elec;Shot.pos],'rows');
[tf,idc]=ismember(N.elec,pos,'rows');
[tf,itt]=ismember(Shot.pos,pos,'rows');
N.elec=pos;
N.a=idc(N.a);N.m=idc(N.m);
fi=find(N.b);N.b(fi)=idc(N.b(fi));
fi=find(N.n);N.n(fi)=idc(N.n(fi));
% showdata2d(N);
Shot.pos=pos;
for i=1:length(Shot.ns), Shot.ns{i}=itt(Shot.ns{i});
    Shot.nx{i}=itt(Shot.nx{i}); end
plotshot(Shot)
[pp,ff,ee]=fileparts(ttfile);
savesgtfile(fullfile(pp,['j' lower(ff) ee]),Shot);
[pp,ff,ee]=fileparts(dcfile);
saveinv2dfile(fullfile(pp,['j' lower(ff) ee]),N);

