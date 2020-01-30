function z=log_odd(r);

% Z=log_odd(R)
%
% this function performs the log odd transform of a population
% in the limit [0;1] into a near normal population. 
% 
% see inverse_log_odd.m for the r to z transform
%
% Vincent MORON
% feb. 2005

z=log(r./(1-r));

