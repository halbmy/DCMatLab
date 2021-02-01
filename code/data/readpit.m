function Shot = readpit(filename)

% READPIT - read pit (picus) file
% Shot = readpit(filename)

Shot=[];Shot.nx={[]};Shot.ns={};Shot.t=[];Shot.tt={[]};Shot.sd=[];
% filename='picus\02 Eucalyptus 180cm.pit';
fid=fopen(filename,'r');
try,
    zeile=fgetl(fid);
    while(~isequal(zeile,'[Main]')), zeile=fgetl(fid); end
    % BPoints are Topo points (up to now they equal MPoints)
    while(~isequal(zeile,'[BPoints]')), % measuring points
        zeile=fgetl(fid);
        [id,vec]=zeil2vec(zeile);
        if isequal(id,'Sensoranzahl'), nsen=vec(1); end
    end
    while(~isequal(zeile,'[MPoints]')), zeile=fgetl(fid); end
    el=zeros(nsen,2);
    for i=1:nsen, %search for measuring points 'i=x/y' in cm
        zeile=fgetl(fid);
        [id,vec]=zeil2vec(zeile);
        el(i,1:2)=vec(1:2);
    end
    Shot.pos=el/100; % in cm
    while(isempty(strfind(zeile,'[oLink'))), 
        zeile=fgetl(fid); end %find first shot
    for i=1:nsen,
        Shot.ns{i}=i;Shot.nx{i}=[];Shot.tt{i}=[];
        zeile=fgetl(fid);
        while(isempty(strfind(zeile,'[o'))), 
            [id,vec]=zeil2vec(zeile);
            vec=vec(find(vec));
            if ~isempty(id),
                j=str2num(id);
                if((~isempty(j))&&(j>0)&&(i~=j)&&(~isempty(vec))), % Shot and receiver different
                    mt=median(vec)/1000; %µs % or mean?
                    st=std(vec)/1000;
                    Shot.nx{i}(end+1)=j;
                    Shot.t(end+1)=mt/1000;
                    Shot.sd(end+1)=st/1000; %/mt
                    Shot.tt{i}(end+1)=mt;
                end
            end
            zeile=fgetl(fid);if isequal(zeile,-1), break; end
        end
    end
    fclose(fid);
catch
    fclose(fid);
    display(lasterr);
end
Shot.t=Shot.t(:);Shot.sd=Shot.sd(:);
Shot.va=[];
for i=1:length(Shot.ns),  for j=1:length(Shot.nx{i}),
      Shot.va(end+1)=sqrt(sum(Shot.pos([Shot.ns{i}(1) Shot.nx{i}(j)]).^2));
    end; end
Shot.va=Shot.va(:)./Shot.t; % apparent velocity