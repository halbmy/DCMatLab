function [model]=invertData( R )

%Ne=readinv2dfile('geomCirc15.elecs');
%Ne=readinv2dfile('geomEllip15_10.elecs');
Ne=readinv2dfile('geomHuman.elecs');
N=readinv2dfile('newDD16.dat');
N.elec=Ne.elec;
am = N.a * 100 + N.m;
ma = N.m * 100 + N.a;
fi1 = find( N.a < N.m );
fi0 = find( N.a > N.m );
[tf, loc] = ismember( am(fi1), ma(fi0) );
fi2 = fi0( loc );
    
system('mkdir tmp');
model=[];
for i=1:size( R, 1 )    
    disp(i)
    R1=R( i, fi1 );
    R2=R( i, fi2 );
    %figure(1);
    %plot(R1,R2,'.');xlabel('normal');ylabel('reverse');
    rez=(R1-R2)./(R1+R2);
    %figure(2);
    %hist(rez*100,30);xlabel('reciprocity/%');ylabel('frequency');
    std_max_rez=[std(rez) max(abs(rez))]*100;

    N.rho=R( i,: );    

    Hin=N; Hin.rho=R1; Hin.k=(fi1);
    Hin.a=N.a(fi1); Hin.b=N.b(fi1); Hin.m=N.m(fi1); Hin.n=N.n(fi1); 

    Rev=N; Rev.rho=R2; Hin.k=(fi2);
    Rev.a=N.a(fi1); Rev.b=N.b(fi1); Rev.m=N.m(fi1); Rev.n=N.n(fi1); 
    
    Med=Hin; Med.rho=(R1+R2)./2; 
    
    if (i==1)
        saveinv2dfile( ['tmp.ohm'], Med); 
        system( 'dc2dtreepre -B tmp.ohm' );
        system('dc2dtreerun tmp.data');
        system('cp `ls -1t tmp/model_* | head -n1` startmodel.vector');
        Nk=readinv2dfile( 'tmp.data' );
    else
        Ncorr = Med;
        Ncorr.konf = Nk.k;
        Ncorr.k = Nk.k;
        Ncorr.r = abs( Med.rho' .* Nk.k );
        Ncorr.rho=0;
        saveinv2dfile( ['tmp.data'] , Ncorr); 
        system('dc2dtreerun -B -r startmodel.vector tmp.data');
    end
    
    system('meshconvert -d2 -A -m -o model -b`ls -1t tmp/model_* | head -n1` -id tmp/meshPara.bms');
    model = cat(3, model, plotInvertData('model.xyr', 'log'));
end
end