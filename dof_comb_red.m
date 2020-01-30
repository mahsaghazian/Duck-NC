function [DOF]=dof_comb_red(x,l,nsim,value,missing,replace)

% [DOF]=dof_comb_red(x,l,nsim,value,missing,replace)
%
% This function computes a combinatorial analysis of the degrees
% of freedom of any input matrix. The total number of combination is based 
% on the number of simulations. This function is useful when n!/(l!(n-l)!).
% is too large. Otherwise, consider dof_comb.m instead.
% 
% Input
% 'x' : matrix of real number (rows describe time and columns describe
% space)
% 'l' : any scalar giving the number of columns choosen each time
% 'nsim' : scalar giving the number of combinations computed.
% 'value' : real number to initiate random seed (if value < 0, it is set
% with the clock of the computer and changes everytime).
% 'missing' : any scalar giving the code of missing entries in x (it could
% be NaN)
% 'replace' : string character to choose to replace missing entries by 
% long-term mean (if replace='yes') or left it as NaN (if replace='no'). In
% the latter case, covariance and dof of 'x' is computed only using available
% entries
%
% Output
% 'DOF' : vector of real number giving the degrees of freedom for every
% combination of 'l' columns of 'x'. The length of the vector is 
% n!/(l!(n-l)!).
%
% Vincent Moron
% Nov. 2005

if(value>0);
    rand('seed',value);
    randn('seed',value);
else
    rand('seed',sum(100*clock));
    randn('seed',sum(100*clock));
end

if missing ~=NaN
  x(find(x==missing))=NaN*ones(size(find(x==missing)));
end

[r,c]=size(x);
for i=1:nsim; N(i,:)=randperm(c); end

if isempty(find(isnan(x)));
    for i=1:nsim;
        DOF(i)=dof(x(:,N(i,1:l)));
    end
    else
    for i=1:nsim;
        DOF(i)=nandof(x(:,N(i,1:l)),replace,missing);
    end
end
    