function showdata(Data)

% SHOWDATA - Show data of any kind or dimension
if isstr(Data), 
    Data=readdatafile(Data); 
end
if isfield(Data,'zweid'), % 3d profile data
    plotprofiles(Data);
else % no 3d profile data
    if ~isfield(Data,'elec')&&isfield(Data,'pos'), Data.elec=Data.pos; end
    if isfield(Data,'x'), % new tokens in
        %version 1: create elec out of it
        
        %version 2: do something real
        %patch2ddata(
        %return;
    end
    if isfield(Data,'elec'),
        %call old showdata
        showdata2d(Data);
        Data.x=Data.elec(:,1);
        if size(Data.elec,2)>1, Data.z=Data.elec(:,end); end
        %patch2ddata(Data)
    else
        display('could not display data!');
    end
end