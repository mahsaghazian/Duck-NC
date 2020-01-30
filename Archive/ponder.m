function [x,pos]=ponder(x,lat);

% [xp]=ponder(x,lat)
%
% This function weights some points according to their
% latitude. The points must be entered
% in column order and the vector describing the latitude
% must be of the same order. The weight for each point is the
% the cosine of latitude.
% 
% Input
% 'x' : matrix of real numbers (row describes time and column describe space)
% containing the data to be weighted.
% 'lat' : vector of real number giving the latitudes of each point of 'x'.
%
% Output
% 'xp' : weighted matrix of real number  
%
% Vincent MORON
% moron@traviata.bo.cnr.it
% (last revision : 02/07/96)

[nl,nc]=size(x);

po=(lat./90).*(pi/2);
pos=cos(po);
for i=1:nc
   x(:,i)=x(:,i).*pos(i);
end


