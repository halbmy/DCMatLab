function rhoa=maketensor(A1,B1,A2,B2,M1,M2,N,RR)

% A1=[3479440.510  5955975.417];
% B1=[3479488.288  5956009.845];
% B2=[3479526.362  5955995.382];
% A2=[3479527.627  5956067.524];
% N=[3479366.392  5955645.778];
% M1=[3479368.670  5955720.747];
% M2=[3479430.598  5955607.125];
% E-Felder E11 E12 E21 E22
% EE=-[0.00003507222685466755 -0.00001853211303934707...
%    -0.00001882070072400068  0.00000071857579408964];

n1=(M1-N);n2=(M2-N);

% E-feld = - delta u / delta x
E11=-RR(1)/norm(n1); % 1. Einspeisung, 2. Messung
E21=-RR(2)/norm(n1);
E12=-RR(3)/norm(n2);
E22=-RR(4)/norm(n2);

% Stromdichtevektoren
[j1x,j1y]=getjvec(A1,N,B1);
[j2x,j2y]=getjvec(A2,N,B2);

% Projektion der E's auf x/y
n1=n1/norm(n1);
n2=n2/norm(n2);

% E1x=(E11*n2(2)-E12*n1(2))/(n1(1)*n2(2)-n2(1)*n1(2));
% E1y=(E11*n2(1)-E12*n1(1))/(n1(2)*n2(1)-n2(2)*n1(1));
% E2x=(E21*n2(2)-E22*n1(2))/(n1(1)*n2(2)-n2(1)*n1(2));
% E2y=(E21*n2(1)-E22*n1(1))/(n1(2)*n2(1)-n2(2)*n1(1));
EE=[n1;n2]\[E11;E12];E1x=EE(1);E1y=EE(2); % geht auch
EE=[n1;n2]\[E21;E22];E2x=EE(1);E2y=EE(2);

jj=j1x*j2y-j1y*j2x;
rhoa=[E1x*j2y-E2x*j1y E2x*j1x-E1x*j2x;...
      E1y*j2y-E2y*j1y E2y*j1x-E1y*j2x]/jj;
% [E1x*j2y E2x*j1y E2x*j1x E1x*j2x;E1y*j2y E2y*j1y E2y*j1x E1y*j2x]/jj
% rhoa=[E1x E2x;E1y E2y]/[j1x j2x;j1y j2y]; % geht auch

if nargout>0, return; end
% Darstellung der Stromdichte für 1. Einspeisung (rot-rot)
plot(A1(1),A1(2),'rx',B1(1),B1(2),'ro',A2(1),A2(2),'gx',B2(1),B2(2),'go',...
    M1(1),M1(2),'mx',N(1),N(2),'mo',M2(1),M2(2),'cx','MarkerSize',10);
axis equal tight
xlabel('rechtswert')
ylabel('hochwert')
grid on
ALL=[A1;A2;B1;B2;M1;M2;N];
dx=10;le=10;
x=min(ALL(:,1))-30*dx:dx:max(ALL(:,1))+30*dx;
y=min(ALL(:,2))-3*dx:dx:max(ALL(:,2))+3*dx;
xlim([min(x) max(x)]);ylim([min(y) max(y)]);
for i=1:length(x),
    for j=1:length(y),
        jvec=getjvec(A1,[x(i) y(j)],B1);
        jabs=sqrt(sum(jvec.^2,2));
        line([x(i) x(i)+jvec(1)/jabs*le],[y(j) y(j)+jvec(2)/jabs*le]);
    end
end
E1a=sqrt(E1x^2+E1y^2);E2a=sqrt(E2x^2+E2y^2);
line([N(1) N(1)+E1x/E1a*le*1],[N(2) N(2)+E1y/E1a*le*1],'color','red')
line([N(1) N(1)+E2x/E2a*le*1],[N(2) N(2)+E2y/E2a*le*1],'color','green')

