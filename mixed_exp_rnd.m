function [X]=mixed_exp_rnd(value,n,alpha,beta1,beta2);

% [X]=mixed_exp_rnd(value,n,alpha,beta1,beta2);
%
% This function computes a random sample of size 'n' from the three parameters
% describing a mixed-exponential distributions
%
% Input
% 'value' : real number to initiate the random sequence (if > 0, the seeds
% is initiated to the value; otherwise, it is changed from the clock)
% 'n' : integer number = length of the time series
% 'alpha' : real number (for choosing one of both means) 
% 'beta1' : real number (= mean of first exponential distribution)
% 'beta2' : real number (= mean of second exponential distribution)
%
% Output
% 'X' : vector of length 'n'
% 
% Vincent MORON
% december 2004

if(value>0);
    rand('seed',value);
    randn('seed',value);
else
    rand('seed',sum(100*clock));
    randn('seed',sum(100*clock));
end

m = [beta1;beta2];           % the two means
t = 1 + (randn(n,1)<alpha);   % randomly pick 1st or 2nd mean
X = m(t) .* -log(rand(n,1)); % generate exponentials with selected means
