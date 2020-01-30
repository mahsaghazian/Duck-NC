function [x_rand]=swg_weibull(value,n,p01,p11,param,nsim);

% [x_rand]=swg_weibull(value,n,p01,p11,param,thres,nsim);
% 
% Stochastic Weather Generator for rainfall (first-order markov chain is 
% used to generate the occurence of rainfall and weibull distribution to 
% simulate rain amounts during wet days). 
%
% Inputs
% 'value' : real number to initiate the random sequence (if > 0, the seeds
% is initiated to the value; otherwise, it is changed from the clock)
% 'n' : integer number = length of simulated time series 
% 'p01' and 'p11' : real number that are respectively the
% probability transition from state 0 to state 1 and the persistance of
% state 1 
% 'param' : vector of two real number giving the scale and shape parameters 
% of the Weibull distribution  
% 'nsim' : integer number = number of simulation. 
%
% Outputs
% 'x_rand' : matrice of n rows and nsim columns giving the simulated
% sequences
%
% Vincent Moron
% November 2004
% revision 1 (nov 2005) : the random is changed into 'seed' method

if(value>0);
    rand('seed',value);
    randn('seed',value);
else
    rand('seed',sum(100*clock));
    randn('seed',sum(100*clock));
end

[x_rand]=markov_chain_1(value,n,p01,p11,nsim);

% simulation of the quantity of rainfall
for i=1:nsim
    a=find(x_rand(:,i)==1);
    ZZ=wblrnd(param(1),param(2),length(a),1);
    x_rand(a,i)=ZZ;
end
        
    

