function Shot1=korrad(Shot)

% KORRAD - corrigate Shot positions to ly on radius if feasible
% Shot_new = korrad(Shot)

mid=mean(Shot.pos);
if max(abs(mid))<0.05, mid=zeros(1,size(Shot.pos,2)); end
rel=Shot.pos-repmat(mid,size(Shot.pos,1),1);
rad=sqrt(sum(rel(:,1).^2+rel(:,2).^2,2));
ang=atan2(rel(:,2),rel(:,1))*180/pi;
Shot1=Shot;
if min(rad)/max(rad)>0.99, % on circle
   ra=rndig(mean(rad),2);
   da=2.5;
   [ang,ii,jj]=unique(round(ang/da)*da);
   Shot1.pos=[ra*cos(ang*pi/180) ra*sin(ang*pi/180)];
   Shot1.pos=Shot1.pos+repmat(mid,size(Shot1.pos,1),1);
   for i=1:length(Shot.ns),
      Shot1.ns{i}=jj(Shot.ns{i});
      Shot1.nx{i}=jj(Shot.nx{i});
   end
end