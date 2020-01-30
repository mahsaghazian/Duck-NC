function [x,y]=persist(p01,p11);

% [x,y]=persist(p01,p11)
%
% this function transforms the p01 and p11 scalar values into x and y.
% x is the probability of a wet day and y is the first order correlation or
% persistence of the daily precipitation occurence. p01 is the probability of
% precipitation, given no precipitation on the previous day and p11 is the 
% probability of precipitation, given precipitation in the previous day.
% The back transform of x,y into p01 and p11 is given by function inverse_persist.m
% p01 and p11 could be computed with function prob_wet_dry.m.
%
% Vincent MORON
% feb 2005

x=p01./(1+p01-p11);
y=p11-p01;