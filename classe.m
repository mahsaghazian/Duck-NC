function [IND,N,FREQ]=classe(X,thres);

% [ind,n,f]=classe(X,thres);
%
% This function classifies the observation in vector "X" according to the
% thresholds defined in 'thres'. The outputs 'ind' is an index with 1 corresponding to
% the first class, 2 corresponding to the second class, etc., 'n' is the number of 
% events in each class and 'f' is the simple frequency in each class. For
% example, if X=[0 0 1 2 5 -2 2 2 55 10 2 2], the command
% [Ind,n,freq]=classe(X,[0 1 2 5 10 20])
% 
% Ind =
%
%     1     1     2     3     4     1     3     3     7     5     3     3
%
% n =
%
%     3     1     5     1     1     0     1
%
% freq =
%
%    0.2500    0.0833    0.4167    0.0833    0.0833         0    0.0833 
%
% Input
% 'X' : vector of real number
% 'thres' : threshold for the classes
% 
% Output
% 'ind' : vector of integer giving the class (from the lowest to the
% highest values) of each unit
% 'n' : vector of integer (of length 'thres'+2) giving the number of unit
% in each class.
% 'f' : vector of real number (of length 'thres'+2) giving the proportion 
% of the number of unit in each class.
%
% Vincent MORON
% March 2005

thres=[-inf thres inf];
l=length(thres)-1;
X=X(:);
X=X(find(~isnan(X)));
[nr]=length(X);

IND=zeros(size(X));
for i=1:l;
    F=find(X > thres(i) & X <= thres(i+1));
    N(i)=length(F);
    FREQ(i)=N(i)./nr;
    IND(F)=i*ones(length(F),1);
end

IND=IND';
