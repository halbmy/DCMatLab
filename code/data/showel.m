function showel(N,start)

% SHOWEL - Show Electrode Configurations
% showel(N[,startwith])

if nargin<2, start=1; end

for i=start:1:length(N.r),
    plot(N.elec(:,1),N.elec(:,2),'b+');
    hold on
    plot(N.elec(N.a(i),1),N.elec(N.a(i),2),'r+');
    plot(N.elec(N.m(i),1),N.elec(N.m(i),2),'g+');
    if N.b(i)>0,
        plot(N.elec(N.b(i),1),N.elec(N.b(i),2),'r+');
    end
    if N.n(i)>0,
        plot(N.elec(N.n(i),1),N.elec(N.n(i),2),'g+');
    end
    hold off
    title(num2str(i))
    pause
end    
