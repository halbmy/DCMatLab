clear all;
clc
%% definitionen

    [fname,pname] = uigetfile('*.txt*');
    filename      = fullfile(pname,fname);

wlength = 70;  % fensterlänge in m
wstep   = 1;   % steplänge in m
ncol    = 2;   % Spalte (Log)
fmin    = .008; % Minimale Frequenz (cycles/m)
fmax    = .9;  % Maximale Frequenz (cycles/m)
zdetail = 200; % tiefe im detail

%% check konsitenz zwischen minimaler frequenz und fensterlaenge
if fmin < (1/(2*wlength))
    disp('minimale frequenz zu klein');
end
%% daten einlesen
A=textread(filename,'','headerlines',2);

%% daten zuordnen
t=A(:,1);    % tiefenvektor
w=A(:,ncol); % messwerte
%plot(t,w)

%% check messwerte
fi=find(w<-999);
if length(fi)>0,
    error('Fehlerwerte!');
end

%% tiefenschrittweite
dt = diff(t);
if max(dt)-min(dt)>0.01,
    warning('Teufenvektor nicht äquidistant!');
end

dt      = median(dt);
ilength = round(wlength/dt); % Fensterlänge in samples
istep   = round(wstep/dt);   % steplänge in samples

%% beginne analyse
w_org=w;
for r=1:2
    if r==1
        % test variance
        % generate test data with variance of the measured data
        w = mean(w_org) + std(w_org)*randn(size(w_org));
    else
        % jetzt messwerte
        w = w_org; % messwerte
    end
    %% frequenz gehalt
    N    = 2^(nextpow2(ilength)+1);
    freq = (1/dt)/N*(1:N/2-1); % frequenz
    wl   = 1./freq;            % wellenlaenge
    
    w_unfiltered = w; % backup
    
    %% bandpass
    samplingFreq = 1/dt;       % 1/tiefenschrittweite
    eckFreq      = [fmin fmax]; % definiert den freqbereich der "drin" bleibt
    
    n     = 2;
    Wn    = eckFreq./samplingFreq;
    ftype = 'bandpass';
    [b,a] = butter(n,Wn,ftype);
    
    w            = filtfilt(b,a,w_unfiltered); % wende filter an
    
    %% berechne FFT
    nwin=floor((length(w)-ilength)/istep);
    F=zeros(length(wl),nwin);
    for i=1:nwin,
        istart = (i-1)*istep+1;
        iend   = istart+ilength-1;
        wi     = detrend(w(istart:iend));
        fi     = fft(wi,N);
        si     = abs(fi(2:N/2));
        F(:,i) = si; % amplitudenspektrum
        FP(:,i)= fi(2:N/2).*conj(fi(2:N/2)); % powerspektrum
    end
    
    if r==1
        Psigma = max(FP'); % signifikanzschwelle aus FFT der randomwerte
    end
end

%% Plot
%figure
%clf;
%fw=find(freq<=fmax & freq>fmin); % nur zwischen fmin  und fmax plotten
fw      = find(freq<=fmax); % nur kleiner fmax plotten
maxamp  = mean(Psigma(fw)); % aus signifikanzschwelle
%maxamp  = 10000000; % manuelle signifikanzschwelle 
nx     = 15; % Wert runter setzen: Label-Frequenz niedriger

clf;
% daten
subplot(3,4,[4 8 12])
    plot(w_unfiltered,t)
    hold on
    %title('Log')
    axis ij; grid on
    xlim(1.1*[min(w_unfiltered) max(w_unfiltered)])
    ax1=gca;
    ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor',[1.0 0.0 0.0],'YColor',[1.0 0.0 0.0]);
    line(w,t,'Marker','none','LineStyle','-','Color',[1.0 0.0 0.0],'Linewidth',2,'Parent',ax2);
    axis ij;
    ylabel('depth (m)')
    linkaxes([ax1 ax2],'y');  
    
% %amplitudenspektrum
% subplot(1,4,[1:3])
%     FW = (F(fw,:));
%     z  = min(t)+ilength*dt/2+(0:nwin-1)*istep*dt;
%     contourf(1:length(fw),z,FW');
%     contourf(freq(fw),z,FW',20);
%     shading flat
%     colorbar
%     axis ij
%     title('Amplitude spectra')

% powerspektrum
subplot(3,4,[1:3 5:7])
    FPW = (FP(fw,:));
    FPW(FPW < maxamp) = maxamp;
    z=min(t)+ilength*dt/2+(0:nwin-1)*istep*dt;
    contourf(freq(fw),z,(FPW'),50);
    %contourf(freq(fw),z,FPW',50);
    shading flat
    %colorbar
    axis ij
    title('Power spectral density vs. depth')
    %xlabel('Frequency (cycles/m)')
    ylabel('depth (m)')
    set(gca,'xlim',[1e-4 fmax])
    set(gca,'xtick',linspace(fmin,fmax,nx));
%     set(gca,'xtick',[.01 .1]);
    
    % achsenbeschriftung auf wellenlaenge
     xt1=rndig(1./get(gca,'xtick'));
     set(gca,'XTickLabel',num2strcell(xt1));
%     for i=1:length(fw1),
%        te=text((i-1)*nx+1,z(1),num2str(rndig(freq(fw1(i)))));
%        set(te,'VerticalAlignment','bottom','HorizontalAlignment','center'); 
%     end




subplot(3,4,[9:11])
    zi=find(z>zdetail,1);
    plot(freq(fw),(FP(fw,zi)))
    hold on
    plot(freq(fw),(Psigma(fw)),'k--')
    grid on
    title(['Power spectral density at ' num2str(z(zi)) 'm'])
    xlabel('wavelength (m)')
    ylabel('')
    set(gca,'xlim',[1e-4 fmax])
     % achsenbeschriftung auf wellenlaenge
     xt2=rndig(1./get(gca,'xtick'));
     set(gca,'XTickLabel',num2strcell(xt2));
    

% %%
% npick=6
% % x-Koordinate eingeben um aus Plot Werte abzulesen (WL und Freq.)
% wl(npick)
% freq(npick)
% %%
% %vline(3,'w','68.3m');
% 
% 
% 
