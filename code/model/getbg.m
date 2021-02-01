function bg=getbg(x,y,M,el)

% GETBG - returns (background) parameter
% bg = getbg(x,y,M,el)

%M=magic(3);x=0:3;y=0:3;el=[2 2];
nel=size(el,1);
l=1;
bg=ones(nel,1);
for l=1:nel, % alle Elektroden
    i1=max(find(x<=el(l,1)));
    i2=min(find(x>=el(l,1)));
    j1=max(find(y<=el(l,2)));
    j2=min(find(y>=el(l,2)));
    if ~isempty(i1*i2*j1*j2),
        ms=M(i2-1:i1,j2-1:j1);
        bg(l)=exp(mean(log(ms(:))));
    end
end