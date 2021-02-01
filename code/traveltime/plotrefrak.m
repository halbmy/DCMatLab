function plotrefrak(xe,ts,t)

% ne=round(sqrt(length(ts)*2))+1;
ne=length(xe);
l=0;
clf;
hold on
for n=1:ne,
    fi=l+(1:ne-n);
    plot(xe(n+1:end),ts(fi),'bx-');
    if nargin>2, plot(xe(n+1:end),t(fi),'ro:'); end
end
hold off