function Savedata =combinesipdata(Data)

% COMBINESIPDATA - Combine SIP data
% Combdata = combinesipdata(Data)

if ~isfield(Data,'sip'), Savedata=Data;return; end

Savedata=Data.sip{1};
if isfield(Data,'ff'), freq=Data.ff; else freq=1:length(Data.sip); end
Savedata.f=ones(length(Savedata.a),1)*freq(1);
for i=2:length(Data.sip),
    sipdata=Data.sip{i};
    olda=length(Savedata.a);
    la=length(sipdata.a);
    fn=fieldnames(Savedata); % may change later
    for j=1:length(fn), % all fields
        if isfield(sipdata,fn{j}),
            ff=getfield(sipdata,fn{j});
            if length(ff)==la,
                fff=getfield(Savedata,fn{j});
                fff(olda+1:olda+la)=ff(:);
                Savedata=setfield(Savedata,fn{j},fff);
            else, % not same length for current frequency
                Savedata=rmfield(Savedata,fn{j});
            end
        elseif isequal(fn{j},'f'), % frequency
                Savedata.f(olda+1:olda+la)=freq(i);
        else, % not available for current frequency
            Savedata=rmfield(Savedata,fn{j});
        end
    end
end
if isfield(Data,'x'), Savedata.x=Data.x; end
if isfield(Data,'y'), Savedata.y=Data.y; 
elseif isfield(Data,'elec'), Savedata.elec=Data.elec; end
