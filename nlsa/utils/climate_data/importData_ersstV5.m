function varargout = importData_ersstV5( DataSpecs )
% IMPORTDATA_ERSSTV5 Read ERSSTv5 monthly data from a single NetCDF file.
%
%* 今回のERSSTv5はv4と同じ単一NetCDF形式なので、v4 importerを使う。
if ~endsWith( DataSpecs.In.file, '.nc' )
    DataSpecs.In.file = [ DataSpecs.In.file '.nc' ];
end

if nargout > 0
    [ varargout{ 1 : nargout } ] = importData_ersstV4( DataSpecs );
else
    importData_ersstV4( DataSpecs )
end
