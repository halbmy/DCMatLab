

fis=find(N.elec(:,2)==0);
fis0=[fis(:);0];
fi0=find(ismember(N.a,fis0)&ismember(N.b,fis0)&ismember(N.m,fis0)&ismember(N.n,fis0));
NN=extractmeasurement(N,fi0);
[mids,seps,ii,kk]=midkonf2d(N);
% feld=feld(fi);
if length(NN.a)>0, showdata2d(NN); end
fiu=find(N.elec(:,2)~=0);
bpos=unique(N.elec(fiu,1));
l=10;
dx=median(diff(unique(N.elec(fis,1))));
if bpos>1, % crosshole
    for i=1:length(bpos),
        fib1=find(N.elec(:,1)==bpos(i));
        dz=median(diff(unique(N.elec(fib,2))));
        for j=i+1:length(bpos),
            fib2=find(N.elec(:,1)==bpos(j));
            dz2=median(diff(unique(N.elec(fib,2))));
            fi=findmess(N,fib1,fib2);
            if isempty(fi),
                if ~isempty(fi0)||(l>10), figure(l); end
                
                l=l+1;
            end
        end
    end
end
for i=1:length(bpos),
    fib=find(N.elec(:,1)==bpos(i));
    fi=findmess(N,fis,fib);
    dz=median(diff(unique(N.elec(fib,2))));
    if ~isempty(fi),
        if ~isempty(fi0)||(l>10), figure(l); end
        for ii=1:length(fi),
           xx= 
        end
        l=l+1;
    end
end

