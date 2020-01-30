function [p01,p11]=inverse_persist(x,y);

% [p01,p11]=inverse_persist(x,y)
%
% this function transforms back the parameters x and y into p01 and p11.
% x is the probability of a wet day and y is the first order correlation or
% persistence of the daily precipitation occurence. p01 is the probability of
% precipitation, given no precipitation on the previous day and p11 is the 
% probability of precipitation, given precipitation in the previous day.
% The transform of p01,p11 into x and y is given by function persist.m
%
% Vincent MORON
% feb 2005

p11=(x.*(1-y))+y;
p01=p11-y;