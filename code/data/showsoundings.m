function showsoundings(N,R)

% SHOWSOUNDINGS - Show data in form of soundings
% showsoundings(N)
% showsoundings(N,model_response)

if ~isfield(N,'eind'), return; end
nn=length(N.eind);
abmi=Inf;abma=0;rmi=Inf;rma=0;
for i=1:nn,
    ab2=N.eind{i}(:,1);
    rhoa=N.eind{i}(:,2);
    abmi=min(abmi,min(ab2));
    abma=max(abma,max(ab2));
    rmi=min(rmi,min(rhoa));
    rma=max(rma,max(rhoa));
end
for i=1:nn,
    subplot(1,nn,i);
    ab2=N.eind{i}(:,1);
    rhoa=N.eind{i}(:,2);
    loglog(rhoa,ab2,'bo:');
    if nargin>1,
        hold on
        loglog(R(N.nr{i}),ab2,'rx:');
        hold off
    end
    set(gca,'Ydir','reverse');
    xlim([rmi rma]);
    ylim([abmi abma]); 
    set(gca,'XTick',[rmi rma],'YTick',[abmi abma],...
        'XTickLabel',{num2str(round(rmi)),num2str(round(rma))},'YTickLabel',{num2str(abmi),num2str(abma)});    
    xlabel('\rho in \Omegam');
    if i==1, ylabel('AB/2 in m'); end
    grid on
    title(sprintf('x=%d',round(N.pos(i))));
end