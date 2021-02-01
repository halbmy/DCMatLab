function showgeolore(Geolore,col)

if nargin<2, col='b'; end
if nargin<3, dt=mod(datenum('0:10',15),1); end
ii=ceil(min(Geolore.gpstime)/dt):floor(max(Geolore.gpstime)/dt);
tt=ii*dt;
lont=interp1(Geolore.gpstime,Geolore.lon,tt);
latt=interp1(Geolore.gpstime,Geolore.lat,tt);
plot(Geolore.lon,Geolore.lat,col,lont,latt,[col '.']);
for i=1:length(tt), 
    
    text(lont(i),latt(i),[' ' datestr(tt(i),15)],'Color',col); 
end
axis equal tight;grid on;
xlabel('x in m');ylabel('y in m');
