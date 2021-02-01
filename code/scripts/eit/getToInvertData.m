function [ data ] = getToInvertData( opt )
    
    if ( nargin == 0 )
        filename ='k05-s001c.get';
    else
        filename = opt;
    end

    curr = 4.7;  % Strom in mA
    gain_digit = 500;
    gainlin = 10^( (40*((gain_digit*5*1.2)/(4096*5.9)-.25)+30)/20);
    ifaktor = gainlin * curr / 1000.0;
    
    fid = fopen( filename, 'r' );
    get = fread( fid, [256 inf], '4*float' );
	get = get';
	fclose( fid );
               
    % alle 208 Messdaten:
           
    allRho = get( :, 1:208 )./ifaktor;

    data = allRho;
end

