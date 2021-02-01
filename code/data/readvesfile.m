function [out]=readvesfile(vesfile)

% READVESFILE - read vertical sounding file
% out = readvesfile(filename)
% filename should contain ab/2 mn/2 rhoa OR
% headerline defining ab/2,mn/2,rhoa or rho 
%   by keywords ab* mn* rhoa* app* rho*
% followed by values

% vesfile='d:\Guenther.T\3d\quorn\dc\q60e60s.ves';
fid=fopen(vesfile,'r');
zeile=fgetl(fid);fclose(fid);
i=1;
while(i<length(zeile)), 
    if isequal(zeile(i:i+1),[32 32]), zeile(i+1)=''; 
    else i=i+1; 
    end
end    
ft=find((zeile==9)|(zeile==32));
iab=0;imn=0;ir=0;irhoa=0;
if any(ft),
    ft(end+1)=length(zeile)+1;
    if ft(1)>1, ft=[0 ft]; end
    for i=1:length(ft)-1,
       ss=lower(zeile(ft(i)+1:ft(i+1)-1));
       if (length(ss)>1)&&isequal(ss(1:2),'ab'), iab=i; end
       if (length(ss)>1)&&isequal(ss(1:2),'mn'), imn=i; end
       if (length(ss)>2)&&isequal(ss(1:3),'app'), irhoa=i; end
       if (length(ss)>2)&&isequal(ss(1:3),'rho'), 
           if (length(ss)>3)&&isequal(ss(4),'a'), irhoa=i; else ir=i; end
       end
       if (length(ss)>=10)&&(isequal(ss(1:10),'resistance')), ir=i; end
    end
end
[data]=textread(vesfile,'','commentstyle','shell','headerlines',double(iab>0));
if iab==0, iab=1;imn=2;irhoa=3; end % not given
ab=data(:,iab);mn=data(:,imn);
if irhoa, 
    rhoa=data(:,irhoa); 
elseif ir,
    k=(ab.^2-mn.^2)./min(ab,mn)*pi/2;rhoa=data(:,ir).*k;
else
    rhoa=ones(size(data,1),1); 
end
out=[ab mn rhoa];
