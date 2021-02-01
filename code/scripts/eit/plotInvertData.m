function [ui, xs, ys] = plotInvertData( filename, opt )

plotstyle='lin';
if nargin==2
   if ( opt == 'log' )
      plotstyle = 'log';
   end
end

data=load(filename);
xs=linspace( min( data(:,1) ), max( data(:,1) ), 64 );
ys=linspace( min( data(:,2) ), max( data(:,2) ), 64 );
[X,Y]=meshgrid(xs,ys);

cmin = 1; cmax = 30;
if ( plotstyle=='log' )
ui=griddata(data(:,1),data(:,2),log10(data(:,3)), X,Y, 'cubic');
cmin = log10(cmin); cmax = log10(cmax);
else 
ui=griddata( data(:,1), data(:,2), data(:,3), X, Y);
end

imagesc(xs,ys,ui)
set(gca,'ydir','normal')
set(gca,'CLIM',[ cmin, cmax ])
axis equal tight
colorbar;

end