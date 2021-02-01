X=0;Y=0;Z=0; % measuring point -> matrices

%% coordinate transformation (x,y,z)->(u,v,w)->(p,q,r)

% equations (21)
g1=(p2-p1)/(q2+q1);
g2=(p2-p3)/(q2+q1);
h1=(p1*q2-p2*q1)/(q2-q1);
h2=(p3*q2-p2*q1)/(q2-q1);
R1=sqrt(p1.^2+q1.^2+r1.^2);
R2=sqrt(p2.^2+q2.^2+r2.^2);
R3=sqrt(p3.^2+q3.^2+r3.^2);
% preparation for eq. (33)
alpha1=atan(g1);alpha2=atan(g2);
sa1=sin(alpha1);ca1=cos(alpha1); %better without sin/cos?
sa2=sin(alpha2);ca2=cos(alpha2);
% sin= +/- tan/sqrt(1+tan^2) cos= +/- 1/sqrt(1+tan^2)
sa1p1=1+sa1+1;sa2p1=1+sa2+1; % spare brackets
pR1=p1+R1;pR2=p2+R2;pR3=p3+R3; % spare brackets
% equation (33)
I1=h1*ca1*ln((p1*sa1+q1*ca1+R1)/(p2*sa1+q2*ca1+R2))...
    -h2*ca2*ln((p3*sa3+q1*ca2+R3)/(p2*sa2+q2*ca2+R2))+q1*ln(pR1/pR3)...
    +2*r1*( atan((q1*ca1+sa1p1*pR1)/r1/ca1)
           -atan((q2*ca1+sa1p1*pR2)/r1/ca1)...
           -atan((q1*ca2+sa2p1*pR3)/r1/ca2)...
           +atan((q2*ca2+sa1p1*pR2)/r1/ca2) );