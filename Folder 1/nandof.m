function [y] = nandof(x,replace,missing)

% DOF : Estimation of spatial degrees of freedom
% 
% [dofs]=dof(X,replace,missing);
%
% this function computes the spatial degrees of
% freedom (dofs) or embedding dimension 
% in one spatio-temporal field (X) with N time
% intervals and M variable in space. The data must
% be close to a normal distribution. dofs characterize
% linearly independent spatial properties in a
% phase space. The dof estimate may help a a rule of thumb for the
% acurrate number PC-EOF that could be used to
% sufficiently describe the variability of the system
% [see Fraedrich et al. J. Clim., 8, 1995 (361-369)]
%
% Input
% X : matrix of real numbers
% missing : scalar (or NaN) to define missing values
% replace : string character to replace missing entries by the long-term
% mean (if replace = 'yes') and left missing entries as NaN (if replace =
% 'no').
%
% Output 
% dofs : scalar giving the number of dof. The dof are computed without the
% missing values
%
% Vincent MORON
% moron@traviata.bo.cnr.it
% last revision (25/7/96);

x(find(x==missing))=NaN*ones(size(find(x==missing)));

[n,m]=size(x);
if strcmp(replace,'no');
    x=nanstan(x,'s',missing,'no');
    [c,d]=eig(nancov(x,missing));
    d=diag(d);
    y=(m.^2)/(sum(d.^2));
else
    x=nanstan(x,'s',missing,'yes');
    [c,d]=eig(cov(x));
    d=diag(d);
    y=(m.^2)/(sum(d.^2));
end
