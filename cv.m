function [PP]=cv(X,DIMS);

%  [PP]=cv(X,DIMS);
%
% this function computes cross-validation periods from
% a nb of time increment. The duration of each subperiod of training
% is given in DIMS. 
% 
% Input
% 'X' : integer giving the total length of the time series 
% 'DIMS' : integer giving the length of training period
%
% Output
% 'PP' : matrix of 0 and 1 (X rows and columns= number of training periods) 
% For each column, "1" defines the verification period and "0" the training
% period
%
% Vincent MORON
% Jan 2001

NB=ceil(X/DIMS);

PP=zeros(X,NB);

for i=1:NB,
	PP(((i-1)*DIMS)+1:i*DIMS,i)=ones(DIMS,1);
end;

PP=PP(1:X,:);