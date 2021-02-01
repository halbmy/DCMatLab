function [mid,ang,ori]=getdatpar(N)

% GETDATPAR - Get (circle) data parameters 
%             midpoint, angle and orientation

mid=round(mean(N.elec)*10000)/10000;
nel=size(N.elec,1);
di=N.elec(:,1:2)-repmat(mid(1:2),nel,1);
an=atan2(di(:,2),di(:,1))*nel/pi/2;
fi=find(an<0);an(fi)=an(fi)+nel;
ang=an(1);
ori=an(2)-an(1);