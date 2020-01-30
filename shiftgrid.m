function X=shiftgrid(Y,lon,lat,choice);

% X=shiftgrid(Y,lon,lat,choice);
%
% This function shift a regular grid with one block in eastern longitudes
% and one in western longitudes in a grid from west to east around 0 degre
% of longitude.
% 
% Input
% 'Y' : a matrix of real number containing time in rows and space in
% columns 
% 'lon' : a vector of real number giving the longitudes of the points in
% west-east format (for example [-30:5:50] for 30W-50E with a step of 5
% degres.
% 'lat' : a vector of real number giving the latitudes of the points.
% 'choice' : character string ='lon-lat' is the first columns of 'X'
% describes all latitudes of the first longitudes and ='lat-lon' if the 
% first columns of 'X' describes all longitudes of the first latitude.
%
% Output
% 'X' : a matrix of real number containing the shifted matrix with rwos
% describing time and columns describing space.
%
% Vincent Moron
% Nov. 2005

[nr,nc]=size(Y);
nlon=length(lon);
nlat=length(lat);
dec=length(find(lon < 0))+1;

for i=1:nr;
    if strcmp(choice,'lon-lat');
        P=reshape(Y(i,:),nlat,nlon);
        P=decalage(P,dec);
        X(i,:)=P(:)';
    elseif strcmp(choice,'lat-lon');
        P=reshape(Y(i,:),nlon,nlat)';
        P=(decalage(P,dec))';
        X(i,:)=P(:)';
    end
end