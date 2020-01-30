function [indice]=choose_area_gridded(lon,lat,flon,flat,choice);

% [indice]=choose_area_gridded(lon,lat,flon,flat,choice);
%
% this function finds indice corresponding to a window in a sample of
% longitudes, latitudes coordinates. It serves to extract a spatial window
% from large gridded dataset (i.e. indice are the columns of the large-scale
% matrix).
%
% Input
% 'lon' : vector of real numbers describing the longitudes of grid-points
% 'lat' : vector of real numbers describing the latitudes of grid-points
% 'flon' : vector of two real numbers giving the limits of the window in
% longitudes
% 'flat' : vector of tow real numbers giving the limits of the window in
% latitudes
% 'choice' : character string ='lon-lat' is the first columns of the
% large-scale matrix describes all latitudes of the first longitudes 
% and ='lat-lon' if the first columns of 'x' describes all longitudes of 
% the first latitude.
%
% Output
% indice : vector giving the column number of the spatial window in the
% large-scale matrix.
%
% Vincent Moron
% Nov 2005

if min(flon) < 0
    lon(find(lon > 180))=lon(find(lon > 180)) - 360;
elseif max(flon) > 180
    lon(find(lon < 0))=lon(find(lon < 0)) + 360;
end

lon=lon(:);
lat=lat(:);
[lon,lat]=meshgrid(lon,lat);
if strcmp(choice,'lon-lat');
    lon=lon(:); lat=lat(:);
elseif strcmp(choice,'lat-lon');
    lon=lon'; lon=lon(:); lat=lat'; lat=lat(:);
end

indice=find((lon >= flon(1) & lon <= flon(2)) & (lat >= flat(1) & lat <= flat(2)));

   