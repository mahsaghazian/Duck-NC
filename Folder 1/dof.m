function [y] = dof(x)

% [dofs]=dof(X);
%
% This function computes the spatial degrees offreedom (dofs) or embedding 
% dimension of one spatio-temporal field (X) with N time intervals in rows
% and M variable in columns. The data must be close to a normal distribution. 
% dofs characterize linearly independent spatial properties in a
% phase space. This function does not allow for missing values. Use
% nandof.m when missing values are present. The dof estimate may help a 
% rule of thumb for the acurrate number PC-EOF that could be used to
% sufficiently describe the variability of the system
% [see Fraedrich et al. J. Clim., 8, 1995 (361-369)]
%
% Input 
% x : matrix of real number (rows should describe time and columns space). 
% No missing value are allowed.  
%
% Output
% dofs : real number giving the degree of freedom
%
% Vincent MORON
% moron@traviata.bo.cnr.it
% last revision (25/7/96);

[n,m]=size(x);
x=stan(x,'s');
[c,d]=eig(cov(x));
d=diag(d);
y=(m.^2)/(sum(d.^2))
