function [x]=replace_missing(x,missing);

% [x]=replace_missing(y,missing);
%
% this function replaces missing entries by the long-term mean of each
% column of 'x'.
% 
% Input
% 'y' : matrix of real number
% 'missing' : scalar used as a code for missing entries (it could be NaN)
%
% Output
% 'x' : matrix of real number where missing entries has been repolaced by
% the long-term mean (it is done independently on each column).
%
% Vincent Moron
% Nov 2005

if missing ~=NaN
  x(find(x==missing))=NaN*ones(size(find(x==missing)));
end

[nr,nc]=size(x);
xm=nanmean(x);

for i=1:nc;
    a=find(isnan(x(:,i)));
    x(a,i)=xm(i)*ones(length(a),1);
end
