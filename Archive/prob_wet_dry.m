function [p00,p01,p10,p11]=prob_wet_dry(sequence,lseason,threshold)

% [p00,p01,p10,p11]=prob_wet_dry(sequence,lseason,threshold)
%
% this function computes p00 (probability of no precipitation, given no
% precipitation the previous day), p01 (probability of precipitation, given
% no precipitation the previous day), p10 (probability of no precipitation,
% given precipitaion the previous day) and p11 (probability of
% precipitation, given precipitation the previous day) on the sequence of
% rainy days in 'sequence' matrix (if it is a matrix, the probabilities
% are computed columnwise) using the 'threshold' scalar to define wet and
% dry days. sequence > threshold are wet days.
%
% Input 
% 'sequence' : matrix of real numbers with row describing time (typically a
% season) and column describing space (station or grid-point)
% 'lseason' ; scalar giving the length of season
% 'threshold' : real number to define wet days
%
% Output
% 'p00' : matrix of persistance from dry to dry state (rows describe years)
% 'p01' : matrix of persistance from dry to wet state
% 'p10' : matrix of persistance from wet to dry state
% 'p11' : matrix of persistance from wet to wet state
%
% Vincent MORON
% feb 2005

[nt,nv]=size(sequence);
nyear=nt/lseason;

for j=1:nyear;
    seq=sequence(((j-1)*lseason)+1:lseason*j,:);
    for i=1:nv,
        s0=seq(1:lseason-1,i);
        s1=seq(2:lseason,i);
        ndry=length(find(s0 <= threshold));
        nwet=length(find(s0 > threshold));
        p00(j,i)=length(find(s0 <= threshold & s1 <= threshold))/ndry;
        p01(j,i)=length(find(s0 <= threshold & s1 > threshold))/ndry;
        p10(j,i)=length(find(s0 > threshold & s1 <= threshold))/nwet;
        p11(j,i)=length(find(s0 > threshold & s1 > threshold))/nwet;
    end
end