function [data, valid] = readSPL( filename, version )

% readSPL - Read resecs binary traces
% [data,valid] = readSPL(filename[,version])

if ( nargin == 1 )
    version = 'neu';
end

valid = 0;
fid=fopen( filename,'r' );
if ( fid < 0 )
  %  [ 'File not found: ' filename ] 
    data = 0;
    return ;
    error( [ 'File not found: ' filename ] ); 
end

valid = 1;
if ( version == 'alt')
    header=fread(fid, 8*16+10, 'char'); 
    data.onTime = fread( fid, 1, 'int32' ); % x008a
    data.delay_ui= fread( fid, 1, 'int32' );
    data.delay_m= fread( fid, 1, 'int32' );
    data.offTime = fread( fid, 1, 'int32' );
    data.duration_m = fread( fid, 1, 'int32' );
    data.nCycles = fread( fid, 1, 'int16' );
    dummy = fread( fid, 4, 'char' );
    allData = fread( fid, 'float' );
else
    header=fread(fid, 7*16+2, 'char'); 
    data.onTime = fread( fid, 1, 'int32' ); % x0072
    data.delay_ui= fread( fid, 1, 'int32' ); 
    data.delay_m= fread( fid, 1, 'int32' );
    data.offTime = fread( fid, 1, 'int32' );
    data.duration_m = fread( fid, 1, 'int32' );
    data.nCycles = fread( fid, 1, 'int16' );
    dummy=fread(fid, 10, 'char');
    allData=fread( fid, 'float' ); %% x008e
end

nData = (data.onTime + data.offTime) * 2 * data.nCycles;
    
fclose(fid);

data.nTraces = floor( size(allData,1) / nData );

data.d=[];
for i=1:data.nTraces
    data.d=[data.d allData( ((i-1)*nData+1):(i*nData), 1) ];
end