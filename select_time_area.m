function [X]=select_time_area(Y,indicet,indices);

% [X]=select_time_area(Y,indicet,indices);
%
% this function extract a spatial and temporal windows in a regular grid
% with the spatial indices in indices and the temporal indices in indicet.
%
% Input
% 'Y' : a matrix of real number containing the whole data with rows
% describing time and columns describing space.
% 'indicet' : a vector of integers giving the temporal indices
% 'indices' : a vectort of integers giving the spatial indices
%
% Output
% 'X' : a matrix of real number containing the selected data with rows
% describing time and columns describing space.
%
% Vincent Moron
% Nov 2005

X=Y(indicet,indices);