function [F]=filtrage(x,choice,N,ord);

% Y = filtrage(X,choice,N,ord)
%
% This function filters the columns of an input matrix with a recursive
% Butterworth filter.
%
% Input
% X : matrix of real numbers (with time in rows and space in columns). The
% columns of 'x' are independetly filtered.
% 'N' : odd scalar giving the order of the filter 
% 'order' : scalar or vector giving the frequency to be filtered in period 
% in unit time (if data are daily data and order = 30, the cut-off is fixed
% at 1/30 cycle-per-day.
% 'choice' : character string to design the filter ('low', 'high', 
% 'bandpass' or 'stop' for low-pass, high-pass, band-pass (frequency between
% limits defined in 'order' are kept), and band-pass (frequency outside
% limits defined in 'order' are kept) filters). if choice = 'bandpass' or
% 'stop', order should be a two-member vector.
% 
% Output
% Y : filtered matrix
%
% created in May 1996 and modified in May 2004
% Vincent Moron
% moron@cerege.fr

[n,nc]=size(x);
ns=round(n/2);
sim1=flipud(x(2:ns+1,:));
sim2=flipud(x(n-ns:n-1,:));
x=[sim1;x;sim2];

if strcmp(choice,'high') | strcmp(choice,'low');
    ord=ord(1);
    order=2./ord;
    disp(['cut-off is set at ',num2str(ord),' time units']);
elseif strcmp(choice,'stop') | strcmp(choice,'bandpass');
    order=2./ord;
    disp(['cut-off is set at ',num2str(ord(1)),' and ',num2str(ord(2)),' time units']);
end

[b,a]=butter(N,order,choice);

for i=1:nc,
    bf(:,i)=filtfilt(b,a,x(:,i));
end;

F=bf(ns+1:n+ns,:);

