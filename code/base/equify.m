function field=equify(field,level,n)

% EQUIFY - restrict field to discrete values
% newfield = equify(field);
% newfield = equify(field,level); % default=1
% newfield = equify(field,level,number); % default=25

if nargin<2, level=1; end
if nargin<3, n=25; end
anz=prod(size(field));

lfield=log(field);

if level==1,
    n=fix(anz/10);
    level=0.01;
    [N,C]=hist(lfield(:),n);
    nada=N>(max(N)*level);
    nlev=fix(length(find(diff(nada)))/2);
    while nlev<anz/50,
        n=fix(n*1.1);
        [N,C]=hist(lfield(:),n);
        nada=N>(max(N)*level);
        nlev=fix(length(find(diff(nada)))/2);
    end
else
    [N,C]=hist(lfield(:),n);
    nada=N>(max(N)*level);
    fprintf('%d levels found!',fix(length(find(diff(nada)))/2));
end

io=1;co=C(1)-(C(2)-C(1))/2;
while 1,
    in=min(find(nada==0));
    if isempty(in), break; end
    eq=max(C(io:in));
    nada(io:in-1)=-1;
    io=min(find(nada==1));
    if isempty(io), return; end
    nada(in:io-1)=-1;
    cn=(C(io-1)+C(in))/2;
    fprintf('Level %.1f-%.1f\n = %.1f\n',exp(co),exp(cn),exp(eq));
    field(find((lfield>=co)&(lfield<=cn)))=exp(eq);
    co=cn;
end
in=length(N);
eq=max(C(io:in));
cn=C(end)+(C(end)-C(end-1))/2;
field(find((lfield>=co)&(lfield<=cn)))=exp(eq);
