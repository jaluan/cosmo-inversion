%
% sf=scalefacs36(sp36,scaling_model)
%
function sf=scalefacs36(sp36,scaling_model)
%
% Check if scaling_model was specified, otherwise set it to default
%
if (~exist('scaling_model','var')), scaling_model = 'all'; end

%
% Setup the scale factors.
%
sf.ST=sp36.ST;
sf.SLth=1;
sf.SLeth=1;
sf.latitude=sp36.latitude;
sf.longitude=sp36.longitude;
sf.elevation=sp36.elevation;
sf.P=sp36.P;

%
% Now that we have pressure, latitude and longitude, we can compute
% the time dependent scaling factors.
%

load pmag_consts
tdsfsample.lat=sf.latitude;
tdsfsample.long=sf.longitude;
tdsfsample.pressure=sf.P;
tdsfsample.elevation=sf.elevation;
tdsfsample.scaling=scaling_model;
sf.tdsf=get_tdsf(tdsfsample,pmag_consts);
