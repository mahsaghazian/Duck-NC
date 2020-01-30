function [x]=nanstan(y,opt,missing,replace);

% [x]=nanstan(y,opt,missing,replace);
%
% This function standardizes the columns of an input matrix to zero mean
% and eventually unit variance. The missing values coded as NaN or any
% missing code are filled by the long-term mean.
% 
% Input
% 'y' : matrix of real number to be standardized
% 'opt : string character ='m' (output data are normalized to zero mean)
% and ='s' (output data are standardized to zero mean and unit variance). 
% 'missing' : scalar defining the missing value (which can be coded either 
% as a number or as NaN).
% 'replace' : string character to replace the missing value by long-term
% mean, that is zero, (if replace='yes') or left as 'NaN (if replace='no');
%
% Output
% 'x' : matrix of standardized data. The missing data are filled with the
% long-term mean or left as NaN
%
% Vincent Moron
% March 1996

y(find(y==missing))=NaN*ones(size(find(y==missing)));

[nl,nc]=size(y);

ym=nanmean(y);
ys=nanstd(y);
ym=ones(nl,1)*ym;
ys=ones(nl,1)*ys;
x=stan(y,opt);

if strcmp(replace,'yes');
   for i=1:nc,
       a=find(isnan(y(:,i)));
       x(a,i)=zeros(length(a),1);
   end
end
  
    
