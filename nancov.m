function [c]=nancov(x,missing)

% [c]=nancov(x,missing)
%
% This function computes the covariance matrix of any input matrice
% containing missing entries. For each pair of columns, the covariance is
% computed only on available observations on both columns.
%
% Input
% 'x' : matrix of real number 
% 'missing' : any real scalar (or NaN)
%
% Output
% 'c' : covariance matrrix of x
%
% Vincent Moron
% nov. 2005

if missing ~= NaN
    x(find(x==missing))=NaN*ones(size(find(x==missing)));
end

[nt,nc]=size(x);
for i=1:nc;
    for j=1:nc;
        a=find(~isnan(x(:,i)) & ~isnan(x(:,j)));
        la=length(a);
        if la > 1;
            c(i,j)=(1/(la-1))*x(a,i)'*x(a,j);
        else
            c(i,j)=NaN;
        end
    end
end
