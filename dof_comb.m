function [DOF]=dof_comb(x,l,missing,replace)

% [DOF]=dof_comb(x,l,missing,replace)
%
% This function computes a combinatorial exhaustive analysis of the degrees
% of freedom of any input matrix. The total number of combination is based 
% on the number of columns (='n') of the input matrix and the dimension 'l' 
% choosen by the user. The number of combination equals n!/(l!(n-l)!).
%
% Input
% 'x' : matrix of real number (rows describe time and columns describe
% space)
% 'l' : any scalar giving the number of columns choosen each time
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
% Nov 2005

[r,c]=size(x);

x(find(x==missing))=NaN*ones(size(find(x==missing)));

if isempty(find(isnan(x)));
    n=combnk([1:c],l);
    for i=1:length(n);
        DOF(i)=dof(x(:,n(i,:)));
    end
else
    n=combnk([1:c],l);
    for i=1:length(n);
        DOF(i)=nandof(x(:,n(i,:)),replace,missing);
    end
end
    
