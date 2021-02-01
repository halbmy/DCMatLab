if 0,
    fid=fopen('core_440m.asc');
    Unit={};Col={};
    zeile=destrip(fgetl(fid));
    while isempty(zeile), zeile=destrip(fgetl(fid)); end
    i=1;[Col{i},zeile]=strtok(zeile);
    while ~isempty(zeile),
        i=i+1;
        [col,zeile]=strtok(zeile);
        if ~isempty(col), Col{i}=col; end
    end
    units=fgetl(fid);
    i=1;[Unit{i},zeile]=strtok(units);
    while ~isempty(units),
        i=i+1;
        [unit,units]=strtok(units);
        if ~isempty(unit), Unit{i}=unit; end
    end
    ss='';for i=1:length(Col), ss=[ss '%f']; end
    A=fscanf(fid,'%f',[length(Col) Inf])';
    % A=mytextscan(fid,ss);
    fclose(fid);
end
A(A<-999)=NaN;
depth=A(:,1);
mindep=122.86;maxdep=224.82;
B=A((depth>=mindep)&(depth<maxdep),:);
mincol=[0 0 0];
maxcol=[5 10000 10000];
figure(1);clf;
ncols=size(B,2)-1;
for ncol=1:ncols,
    logi=B(:,ncol+1);
    logi((logi<mincol(ncol))|(logi>maxcol(ncol)))=NaN;
    B(:,ncol+1)=logi;
    subplot(1,ncols,ncol);plot(logi,B(:,1),'-');axis ij tight;grid on;
    xlabel([strrep(Col{ncol+1},'_',' ') ' in ' Unit{ncol+1}]);ylabel('depth in m');
end
dx=median(diff(B(:,1)));
figure(2);clf;
if diff(minmax(diff(B(:,1))))<0.01,
    for ncol=1:ncols,
        logi=B(:,ncol+1);
        fi1=find(isfinite(logi));    
        fi2=find(isnan(logi));
        logi(fi2)=interp1(B(fi1,1),logi(fi1),logi(fi2),'nearest');
        [rc,ic,ft]=fouriertransform(logi,dx);
        ac=sqrt(rc.^2+ic.^2);
        subplot(ncols,1,ncol);plot(ft,ac);ylim([0 ac(7)]);
        grid on;
        xlabel('wavenumber in 1/m');ylabel('frequency');
        legend(strrep(Col{ncol+1},'_',' '));
        %     loglog(ft,ac)
    end
end
%     ff=fft(logi);
%     plot(ft,abs(ff));ylim([0 abs(ff(5))]);


%     x=min(depth):dx:max(depth);
%     all=harmfit(depth(fi),logi(fi),0,x);
%     ff=fft(all);    

