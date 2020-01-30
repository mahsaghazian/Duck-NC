function [Z]=proj(x,lon,lat,lonout,latout,ch,choice,choice_out,method);

% [Y]=proj(X,lon,lat,lonout,latout,ch,choice,choice_out,method);
%
% This function interpolates spatially the row of an input matrix x (with 
% rows describing time and columns describing space).
% 
% Input
% 'X' : input matrix containing the values to be interpolated. The rows 
% describe time and the columns describe space. THe interpolation is
% performed independently on each row
% 'lon' : vector containing the longitudes of Y
% 'lat' : vector containing the latitudes of Y
% 'lonout' : vector containing the longitudes of the output grid
% 'latout' : vector containing the latitudes of the output grid
% 'ch' : character string (='shift' or 'no' ; if ch='shift', a dataset
% beginning at 0 is shifted by 180 in the output matrix. This is useful
% for interpolating dataset around 0)
% 'choice' : character string ='lon-lat' is the first columns of 'x'
% describes all latitudes of the first longitudes and ='lat-lon' if the 
% first columns of 'x' describes all longitudes of the first latitude.
% 'choice_out' : character string to define the ordering of output matrix.
% ='lon-lat' is the first columns of 'Y' describes all latitudes of the first
% longitudes and ='lat-lon' if the first columns of 'Y' describes all 
% the longitudes of the first latitude.
% 'method' : method to interpolate the data ('linear', 'cubic', 'spline',
% 'nearest') from interp2.m.  
%
% Output 
% 'Y' : interpolated matrix ordered as 'lon-lat' or 'lat-lon'
%
% Vincent Moron
% July 2001
% modified and checked November 2005

[nt,nv]=size(x);
nlon=length(lon);
nlat=length(lat);

if strcmp(ch,'shift');
    lon(find(lon > 180))=lon(find(lon > 180))-360;
    b=find(diff(lon) < 0);
    b=nlon-b;
else
    b=0;
end

lon=decalage(lon(:)',b);
[lon,lat]=meshgrid(lon(:),lat(:));
[lono,lato]=meshgrid(lonout(:),latout(:));

for i=1:nt;
    if strcmp(choice,'lon-lat');
        y=reshape(x(i,:),nlat,nlon);
        if b > 0;
            y=decalage(y,b);
        end
        z=interp2(lon,lat,y,lono,lato,method);
        if strcmp(choice_out,'lon-lat');
            Z(i,:)=z(:)';
        elseif strcmp(choice_out,'lat-lon');
            z=z'; Z(i,:)=z(:)';
        end
    elseif strcmp(choice,'lat-lon');
        y=reshape(x(i,:),nlon,nlat)';
        if b > 0;
            y=decalage(y,b);
        end
        z=interp2(lon,lat,y,lono,lato,method);
        if strcmp(choice_out,'lon-lat');
            Z(i,:)=z(:)';
        elseif strcmp(choice_out,'lat-lon');
            z=z'; Z(i,:)=z(:)';
        end
    end
end
