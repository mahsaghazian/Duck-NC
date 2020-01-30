function [X2]=decalage(X,va);

% [Y]=decalage(X,nc);
%
% This function  shifts the columns of an input matrix 'X' by 'nc' columns.
% 
% Input
% 'X' : matrix of real numbers
% 'nc' : integer giving the number of columns to be shifted. If nc is
% positive, the shift is processed on the right and if nc is negative, the
% shift is processed on the left (for example, if X is a 10 by 5 matrix
% and nc=2, the function will place the last two columns of X in the first
% two columns of Y.
%
% Output 
% 'Y' : shifted matrix
%
% created by JM Gaillard
% March 1993

[li col]=size(X);
if (abs(va)>=col)
  va=va-(fix(va/col)*col)
end
X2=X';
if (va>0)
  for i=1:va;
    X2=[X2(col,:);X2(1:col-1,:)];
  end;
end;
if (va<0)
  for i=va:1:-1
    X2=[X2(2:col,:);X2(1,:)];
  end;
end;
X2=X2';

