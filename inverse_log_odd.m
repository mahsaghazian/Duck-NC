function r=inverse_log_odd(z);

% R=inverse_log_odd(Z)
%
% this function performs the inverse of a log odd transform
% and ouput a transformed vector R in the limit [0;1].
%
% see log_odd.m for the z to r transform
%
% Vincent MORON
% feb. 2005

r=exp(z)./(exp(z)+1);

