function [robs,p]=corr_mc(a,b,nsim,value);

% [robs,p]=corr_mc(a,b,nsim,value);
%
% This function computes the one-sided level of significance for
% correlation between two vectors 'a' and 'b' of the same length. The
% function creates 'nsim pair of correlation between two different set of
% random time series having the same power spectra but random-phases as 'a'
% and 'b' respectively. The level of significance is the proportion of
% random correlations lower (if the observed correlation is negative) or
% higher (if the observed correlation is positive) than the observed
% correlation.
%
% Input
% 'a' : vector of real number giving the first time series
% 'b' : vector of real number giving the second time series
% 'nsim' : integer number giving the number of time series to simulate
% 'value' : real number to initiate the random sequence (if > 0, the seeds
% is initiated to the value; otherwise, it is changed from the clock)
%
% Output
% 'robs' : real number giving the observed correlations between 'a' and 'b'
% 'p' : real number giving the one-sided level of significance of 'robs'.
%
% Vincent Moron
% Sept 2001

[nl]=length(a);
nr=round(sqrt(nsim));
if nargin==3;
    value=-1;
end
as=ebisuzaki(a,nr,value);
bs=ebisuzaki(b,nr,value);
robs=(1/(nl-1))*stan(a,'s')'*stan(b,'s');
rsim=(1/(nl-1))*stan(as,'s')'*stan(bs,'s');
rsim=rsim(:);
if robs < 0
   p=length(find(rsim < robs))/(nr.^2);
else
   p=length(find(rsim > robs))/(nr.^2);
end


   
