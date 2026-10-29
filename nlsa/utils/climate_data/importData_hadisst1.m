function Data = importData_hadisst1( DataSpecs )
% IMPORTDATA_HADISST1 Read HadISST1 monthly SST from a single NetCDF file.
%
%* MATLAB の netcdf.getVar では sst を lon x lat x time として読む。
%* その後の処理と保存形式は ERSST 系 importer と同じにする。

In     = DataSpecs.In;
Out    = DataSpecs.Out;
Time   = DataSpecs.Time;
Domain = DataSpecs.Domain;
Opts   = DataSpecs.Opts;

if ~isfield( Opts, 'ifDetrend' )
    Opts.ifDetrend = false;
end
if ~isfield( Opts, 'ifOutputData' )
    Opts.ifOutputData = true;
end

if Opts.ifCenter && Opts.ifCenterMonth
    error( [ 'Global and monthly climatology removal cannot be ' ...
             'simultaneously selected' ] )
end
ifClim = Opts.ifCenter || Opts.ifCenterMonth;

if Opts.ifCenter
    fldStr = [ Out.fld 'a' ];
else
    fldStr = Out.fld;
end
if Opts.ifCenterMonth
    fldStr = [ Out.fld 'ma' ];
else
    fldStr = Out.fld;
end
if Opts.ifDetrend
    fldStr = [ fldStr 'l' ];
end
if Opts.ifWeight
    fldStr = [ fldStr 'w' ];
end
if Opts.ifAverage
    fldStr = [ fldStr 'av' ];
end
if Opts.ifNormalize
    fldStr = [ fldStr 'n' ];
end
if ifClim
    fldStr = [ fldStr '_' Time.tClim{ 1 } '-' Time.tClim{ 2 } ];
end

dataDir = fullfile( Out.dir, ...
                    fldStr, ...
                    [ sprintf( 'x%i-%i',  Domain.xLim ) ...
                      sprintf( '_y%i-%i', Domain.yLim ) ...
                      '_' Time.tLim{ 1 } '-' Time.tLim{ 2 } ] );
if Opts.ifWrite && ~isfolder( dataDir )
    mkdir( dataDir )
end

limNum = datenum( Time.tLim, Time.tFormat );
startNum = datenum( Time.tStart, Time.tFormat );
nT = months( limNum( 1 ), limNum( 2 ) ) + 1;
idxT0 = months( startNum, limNum( 1 ) );
if ifClim
    climNum = datenum( Time.tClim, Time.tFormat );
    nTClim = months( climNum( 1 ), climNum( 2 ) ) + 1;
    idxTClim0 = months( startNum, climNum( 1 ) );
end

ncId   = netcdf.open( fullfile( In.dir, In.file ) );
idTime = netcdf.inqVarID( ncId, 'time' );
idLon  = netcdf.inqVarID( ncId, 'longitude' );
idLat  = netcdf.inqVarID( ncId, 'latitude' );
idFld  = netcdf.inqVarID( ncId, In.var );

time = netcdf.getVar( ncId, idTime );
nTTot = numel( time );

if idxT0 < 0
    idxT0Read = 0;
    preDeficit = abs( idxT0 );
    disp( sprintf( ...
        'Date range requested preceeds the available date range. Setting the first %i samples to zero.', ...
        preDeficit ) )
else
    idxT0Read = idxT0;
    preDeficit = 0;
end
postDeficit = idxT0Read + nT - preDeficit - nTTot;
postDeficit = max( postDeficit, 0 );
if postDeficit > 0
    disp( sprintf( ...
        'Date range requested exceeds the available date range. Setting the last %i samples to zero.', ...
        postDeficit ) )
end
nTRead = nT - preDeficit - postDeficit;

lon = double( netcdf.getVar( ncId, idLon ) );
lat = double( netcdf.getVar( ncId, idLat ) );
lon = lon( : );
lat = lat( : );
lon( lon < 0 ) = lon( lon < 0 ) + 360;
[ lon, idxSrt ] = sort( lon, 'ascend' );
nX = length( lon );
nY = length( lat );
[ X, Y ] = ndgrid( lon, lat );

fldRead = netcdf.getVar( ncId, idFld, [ 0 0 idxT0Read ], [ nX nY nTRead ] );
fldRead = fldRead( idxSrt, :, : );

isGood = isfinite( fldRead );
for attName = { '_FillValue' 'missing_value' }
    try
        missingValue = netcdf.getAtt( ncId, idFld, attName{ 1 } );
        isGood = isGood & fldRead ~= missingValue;
    catch
    end
end
%* HadISST1 marks 100% sea-ice-covered SST gridboxes as -1000, but this is
%* not the NetCDF _FillValue. Treat those nonphysical SST values as missing.
isGood = isGood & fldRead > -100;
ifXY = X >= Domain.xLim( 1 ) & X <= Domain.xLim( 2 ) ...
     & Y >= Domain.yLim( 1 ) & Y <= Domain.yLim( 2 ) ...
     & all( isGood, 3 );
iXY = find( ifXY( : ) );
nXY = length( iXY );

fldRead = reshape( fldRead, [ nX * nY nTRead ] );
fld = zeros( nXY, nT );
fld( :, 1 + preDeficit : nTRead + preDeficit ) = fldRead( iXY, : );

if Opts.ifWeight
    dLon = abs( [ diff( lon ); lon( 1 ) + 360 - lon( end ) ] ) * pi / 180;
    dLat0 = median( abs( diff( lat ) ) ) * pi / 180;
    dLat = dLat0 * cos( lat * pi / 180 );
    w = dLon * dLat';
    w = w( ifXY );
    w = sqrt( w / sum( w ) * nXY );
    fld = fld .* w;
end

if Opts.ifCenter
    cli = netcdf.getVar( ncId, idFld, [ 0 0 idxTClim0 ], [ nX nY nTClim ] );
    cli = cli( idxSrt, :, : );
    cli = reshape( cli, [ nX * nY nTClim ] );
    cli = cli( iXY, : );
    if Opts.ifWeight
        cli = cli .* w;
    end
    cli = mean( cli, 2 );
    fld( :, 1 + preDeficit : end - postDeficit ) = ...
        fld( :, 1 + preDeficit : end - postDeficit ) - cli;
end

if Opts.ifCenterMonth
    cliData = netcdf.getVar( ncId, idFld, [ 0 0 idxTClim0 ], [ nX nY nTClim ] );
    cliData = cliData( idxSrt, :, : );
    cliData = reshape( cliData, [ nX * nY nTClim ] );
    cliData = cliData( iXY, : );
    if Opts.ifWeight
        cliData = cliData .* w;
    end
    cli = zeros( nXY, 12 );
    for iM = 1 : 12
        cli( :, iM ) = mean( cliData( :, iM : 12 : end ), 2 );
    end
    idxM0 = month( customDatemnth( limNum( 1 ), preDeficit ) );
    for iM = 1 : 12
        idxM = mod( idxM0 + iM - 2, 12 ) + 1;
        fld( :, iM + preDeficit : 12 : end - postDeficit ) = ...
              fld( :, iM + preDeficit : 12 : end - postDeficit ) ...
            - cli( :, idxM );
    end
end

netcdf.close( ncId );

if Opts.ifDetrend
    t = [ 0 : nT - 1 ] / 12;
    beta = zeros( nXY, 2 );
    for j = 1 : nXY
        beta( j, : ) = polyfit( t, fld( j, : ), 1 );
    end
    fld = fld - beta( :, 2 ) - beta( :, 1 ) * t;
end

if Opts.ifAverage
    fld = mean( fld, 1 );
    if ifClim
        cli = mean( cli, 1 );
    end
end

nD = size( fld, 1 );

if Opts.ifNormalize
    l2Norm = norm( cli( : ), 2 ) / sqrt( nTClim );
    fld = fld / l2Norm;
end

gridVarList = { 'lat', 'lon', 'ifXY', 'fldStr', 'nD' };
if Opts.ifWeight
    gridVarList = [ gridVarList 'w' ];
end
if Opts.ifWrite
    save( fullfile( dataDir, 'dataGrid.mat' ), gridVarList{ : }, '-v7.3' )
end

x = double( fld );
varList = { 'x' 'idxT0' };
if ifClim
    varList = [ varList 'cli' 'idxTClim0' 'nTClim' ];
end
if Opts.ifNormalize
    varList = [ varList 'l2Norm' ];
end
if Opts.ifWrite
    save( fullfile( dataDir, 'dataX.mat' ), varList{ : }, '-v7.3' )
end

if nargout > 0
    varList = [ varList gridVarList ];
    if ~Opts.ifOutputData
        varList = varList( 2 : end );
    end
    nVar = numel( varList );
    vars = cell( 1, nVar );
    for iVar = 1 : nVar
        vars{ iVar } = eval( varList{ iVar } );
    end
    Data = cell2struct( vars, varList, 2 );
end
