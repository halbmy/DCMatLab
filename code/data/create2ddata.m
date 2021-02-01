function N=create2ddata(arr,nel,del,eel,seps,bigd)

% CREATE2DDATA - Create 2d surface data set
% Data = create2ddata(arr,nel,del,fel,seps,bigd)
% arr..array(1-Wenner,2-Pole-Pole,3-Dipole-Dipole,4/5-Wenner-beta/gamma,
%            6-Pole-Dipole,7-Schlumberger,8-PD-reverse,9-circulating DD)
% nel..number of electrodes
% del..electrode distance
% fel..first electrode
% sep..separations ([min max])
% big..dipole enlargement ([every shift])

% initialization
if nargin<6, bigd=[0 0]; end
if (nargin<5)||(isempty(seps)), seps=[1 100]; end
if nargin<4, eel=0; end
if nargin<3, del=1; end
if nargin<2, nel=41; end
if nargin<1, arr=1; end
nmin=seps(1);nmax=seps(2);
enl=0;
if min(bigd)>0,
    enl=1;
    every=bigd(1);
    shift=bigd(1);
end
% creation
N=[];
N.elec=(0:nel-1)'*del+eel;N.elec(1,2)=0;
N.a=[];N.b=[];N.m=[];N.n=[]; 
for n=nmin:nmax,
    df=1;ddf=1; % dipole enlargement for high separations
    if enl&ismember(arr,[3 6 7 8 9]), % DD,PD(f/r),CCPPC,SCHLUM
        %             df=fix(n/(every+1))+1;
        %             ddf=fix(n/(shift+1))+1;
        if every>0, df=ceil(n/every); end
        if shift==0, ddf=1; else ddf=ceil(n/shift); end
    end
    first=(1:ddf:nel)';
    abmn=zeros(length(first),4);
    abmn(:,1)=first;
    switch arr,
        case 1, %Wenner
            abmn(:,3)=abmn(:,1)+n;
            abmn(:,4)=abmn(:,1)+2*n;
            abmn(:,2)=abmn(:,1)+3*n;
        case 2, %Pole-Pole
            abmn(:,3)=abmn(:,1)+n*df;
        case {3,9}, %Dipole-Dipole and CC-PP-C
            abmn(:,2)=abmn(:,1)+df;
            abmn(:,3)=abmn(:,2)+n;
            abmn(:,4)=abmn(:,3)+df;
        case 4, %Wenner-beta
            abmn(:,2)=abmn(:,1)+n;
            abmn(:,3)=abmn(:,1)+2*n;
            abmn(:,4)=abmn(:,1)+3*n;
        case 5, %Wenner-gamma
            abmn(:,3)=abmn(:,1)+n;
            abmn(:,2)=abmn(:,1)+2*n;
            abmn(:,4)=abmn(:,1)+3*n;
        case 6, %Pole-Dipole
            abmn(:,3)=abmn(:,1)+n;
            abmn(:,4)=abmn(:,3)+df;
        case 8, %Pole-Dipole reverse
            abmn(:,1)=nel+1-first;
            abmn(:,3)=abmn(:,1)-n;
            abmn(:,4)=abmn(:,3)-df;
            abmn(find(abmn(:,4)<1),4)=-1;
        case 7, %Schlumberger
            abmn(:,3)=abmn(:,1)+n;
            abmn(:,4)=abmn(:,3)+df;
            abmn(:,2)=abmn(:,4)+n;
    end % switch
    fi=find(max(abmn')>nel);
    abmn(fi,:)=[];
    fi=find(min(abmn')<0);
    abmn(fi,:)=[];
    if ~isempty(abmn),
        N.a=[N.a;abmn(:,1)];
        N.b=[N.b;abmn(:,2)];
        N.m=[N.m;abmn(:,3)];
        N.n=[N.n;abmn(:,4)];
    end
end % n loop
if arr==9, % circulating dipole
    nradd=nel-3;
    N.a=[N.a;ones(nradd,1)];
    N.b=[N.b;ones(nradd,1)*nel];
    N.m=[N.m;(1:nradd)'+1];
    N.n=[N.n;(1:nradd)'+2];
end 
N.k=getkonf2d(N);