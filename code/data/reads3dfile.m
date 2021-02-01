function N=readsnd3d(filename)

% READSND3D - Read 3d file by use of soundings
% N = readsnd3d(filename);
% Format:
% sounding_filename_1 x_location y_location angle [flag]
% ...
% sounding files consists of 3 columns for AB/2 MN/2 and apparent resistivity

N=[];
[fpath,fname,fext]=fileparts(filename);
fid=fopen(fullfile(fpath,[fname fext]),'r');
if isequal(fid,-1), error('Could not open file!'); end
try
    zeile=fgetl(fid);
    zeile(zeile>150)='';
    l=0;
    while ischar(zeile),
        l=l+1;
        ishalf=sscanf(zeile,'%*s %*s %*s %*s %d');
        if isempty(ishalf), ishalf=0; else zeile=zeile(1:end-1); end
        vesfile=sscanf(zeile,'%s %*s %*s %*s');
        xyphi=sscanf(zeile,'%*s %f %f %f');
        xpos=xyphi(1);ypos=xyphi(2);phi=xyphi(3);
        fprintf('Sounding %s at x=%.1fm,y=%.1fm,phi=%d\n',vesfile,xpos,ypos,round(phi));
        vesfull=fullfile(fpath,vesfile);
        if ~exist(vesfull,'file'),
            fprintf('VES file %s does not exist!\n',vesfile);
        else
            data=readvesfile(vesfull);
            AA=data(:,1);MM=data(:,2);
            BB=-AA;NN=-MM;
            V2=abmn2n(AA,BB,MM,NN);V=V2;
            if ishalf,
                V.b(:)=0;
            end
            V.elec(:,1)=V2.elec(:,1)*sin(phi*pi/180)+xpos;
            V.elec(:,2)=V2.elec(:,1)*cos(phi*pi/180)+ypos;
            V.elec(:,3)=0;V.elec=round(V.elec*100)/100;
            V.r=data(:,3);%data{3};
            if isempty(N),
                N=V;
            else
                N=combdata3d(N,V);
            end
            eind{l}=data(:,[1 3]);%data{1} data{3}];
            zweid{l}=V2;
            names{l}=strrep(vesfile,'.ves','');
            nr{l}=length(N.r)-length(V.r)+1:length(N.r);
            pos{l}=[xpos ypos phi];
        end
        zeile=fgetl(fid);
    end
    N.zweid=zweid;N.nr=nr;
    N.eind=eind;
    N.names=names;
    N.pos=pos;
    fclose(fid);
catch
    fclose(fid);
    display(lasterr);
    return;
end
if ~isempty(N), N.k=getkonf(N); end
