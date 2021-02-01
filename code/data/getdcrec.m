% function rec = getdcrec(Shot)

% GETDCREC - Get DC reciprocity
% rec = getdcrec(Data)
% Shot..shot structure with pos,ns,nx and t
% rec..relative reciprocity

A=zeros(length(N.a),4);
A(:,1:2)=sort([N.a(:) N.b(:)],2);
A(:,3:4)=sort([N.m(:) N.n(:)],2);
fi=find(N.a>N.m);
du=A(fi,1:2);A(fi,1:2)=A(fi,3:4);A(fi,3:4)=du;
abmn=A(:,1)*999999+A(:,2)*9999+A(:,3)*99+A(:,4);
rec=zeros(l,1);
for i=1:length(N.a),
    fi=find(abmn==abmn(i));
    if any(fi), rec(i)=std(N.r(fi))/mean(N.r(fi)); end
end
