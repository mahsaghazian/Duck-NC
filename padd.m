function [Y]=padd(X,lseason,n);

% [Y]=padd(X,lseason,n);
%
% This function serves to add 'n' data at the beginning and the end of a
% vector/matrix entered in 'X'. 'lseason' is the length of each season of
% 'X' in which rows describing time (daily time scale with seasons beneath
% each others). 'n' is the number of times the first and last observations
% of each season is copied in 'Y';
% if X contains for example 92 days by 39 seasons, the function
% [Y]=padd(X,92,10); will creates a new matrix of 112*39 rows where
% th first day of each season is repeated 10 times at the beginning and the
% last day of each season is repeated 10 times at the end of each season.
%
% created by: Vincent Moron
% August 2005

[nt,nv]=size(X);
nyear=nt/lseason;
lseason2=lseason+(n*2);

for i=1:nyear;
    sample=X(((i-1)*lseason)+1:i*lseason,:);
    Y(((i-1)*lseason2)+1:i*lseason2,:)=[copy(sample(1,:),n);sample;copy(sample(lseason,:),n)];
end

