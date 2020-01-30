function [analog,output]=find_analog(X,library,target,knn,value); 

% [analog,output]=find_analog(X,library,target,knn,value); 
% 
% subroutine for finding the knn analog of a column vector entered
% in 'X' in a library entered in 'library'. The closest analog has more 
% chance to be selected than the second one, which ahve
% more chance to be selected that the third one etc.). The distance is the
% squared Euclidean distance but is should be expanded to other solutions 
% in next versions.
%
% Input
% 'X' : row vector of real number 
% 'library' : matrix  of real number (same number of column as 'X') in
% which analog are searched
% 'target' : matrix of real number (same number of rows as 'library') that
% are associated with 'library'
% 'knn' : scalar integer giving the number of knn selected
% 'value' : real number to initiate the random sequence (if > 0, the seeds
% is initiated to the value; otherwise, it is changed from the clock). By
% default, the latter solution is considered.
% 
% Output 
% 'analog' : vector of integers giving the rank of analogs in 'library' and
% 'target'. 
% 'output' : matrix of real number (= analogs in 'target')
%
% Vincent MORON
% version 1 ; 12/1/2005

if nargin==4;
    value=-1;
end

if(value>0);
    rand('seed',value);
    randn('seed',value);
else
    rand('seed',sum(100*clock));
    randn('seed',sum(100*clock));
end

[nl,nv]=size(library);
nb=floor(sqrt(nl));
p=cumsum((1./[1:nb])./sum(1./[1:nb])); % (cf. Beersma and Buishand, CR, 25, 121-133, 2003)

[nx,vx]=size(X);
[ny,vy]=size(library);
DD=distance_euclid(library,X,'sqeuclid');
[A,B]=sort(DD);
A=A(1:nb); B=B(1:nb);
if ~isempty(knn);
    for k=1:knn;
        uni=rand(1);
        [aa,bb]=sort([uni,p]);
        pos(k)=find(aa==uni);
    end
    analog=B(pos);
    output=target(analog,:);
else
    analog=B(1);
    output=target(analog,:);
end