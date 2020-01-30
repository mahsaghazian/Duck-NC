function [X]=fliplonlat(Y,lon,lat,choicein,choiceout);

% [X]=fliplonlat(Y,choicein,choiceout);
%
% This function serves to flip an input grid from "lon-lat" (where the
% first columns describe all latitudes of the first longitudes) to 
% "lat-lon" (where the first columns describe all longitudes of the first
% latitudes). 
%
% Input
% ''Y' : a matrice of real number containing the data to flip with rows
% describing time and columns describing space
% 'lon' : a vector of real number giving the longitudes of the grid
% 'lat' : a vector of real number giving the latitudes of the grid
% 'choicein' : character string describing the format of the input grid
% 'choiceout' : character string describing the format of the output grid
%
% Vincent Moron
% Nov. 2005

[nr,nc]=size(Y);
nlon=length(lon);
nlat=length(lat);

for i=1:nr;
    if strcmp(choicein,'lon-lat') & strcmp(choiceout,'lat-lon');
        y=(reshape(Y(i,:),nlat,nlon)');
        X(i,:)=y(:)';
    else
        y=(reshape(Y(i,:),nlon,nlat)');
        X(i,:)=y(:)';
    end
end
