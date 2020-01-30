function [d]=distance_euclid(a,b,opt);

% [d]=distance_euclid(a,b,opt);
%
% function to compute the euclidean distance between the rows of an input 
% matrix entered in 'a' and a column-vector entered in 'b'
%
% Input 
% 'a' : matrix of real number
% 'b' : column-vector of real number (same nb of column as in 'a')
% 'opt' : character string to define the ED (='euclid' for normal ED;
% ='sqeuclid' for Squared ED and ='seuclid' for standardized ED.
%
% Ouput
% 'd' : vector of real number giving the ED between the rows of 'a' and
% 'b'.
%
% Vincent Moron
% Nov 2005

[nr,nc]=size(a);
b=copy(b,nr);
v=diag(1./var([a;b(1,:)]));
if strcmp(opt,'euclid'); % Euclidean distance
    if nc > 1
        d=sqrt(sum(((a-b).^2)'));
    else
        d=sqrt(((a-b).^2)');
    end
elseif strcmp(opt,'sqeuclid'); %Squared Euclidean distance
    if nc > 1
        d=(sum(((a-b).^2)'));
    else
        d=(((a-b).^2)');
    end
elseif strcmp(opt,'seuclid'); % Standarddized Euclidean distance
    if nc > 1
        d=sqrt(sum((((a-b).^2)*v)'));
    else
        d=sqrt((((a-b).^2)*v)');
    end
end
