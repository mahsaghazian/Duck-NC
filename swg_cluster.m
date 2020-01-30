function [Y]=swg_cluster(value,seq,nsim,out);

% [Y]=swg_cluster(value,x,nsim,out);
%
% this function creates ranodm time series matching a given probability
% transition.
%
% Input
% 'value' : real number to initiate the random sequence (if > 0, the seeds
% is initiated to the value; otherwise, it is changed from the clock)
% 'seq' : matrice or row-vector of integer number describing for example
% the sequence of weather types in day (in columns) for different years.
% 'nsim' : integer number = number of time series simulated
% 'out' : character string (= name of the output file)
%
% Output
% 'Y' : matrice of n rows (= number of columns in 'seq') and columns (=
% 'nsim). These sequences matches the relative frequency of each integer
% in the input and also the probability transition (order 1) between them.
%
% Vincent MORON
% September 2005

if(value>0);
    rand('seed',value);
    randn('seed',value);
else
    rand('seed',sum(100*clock));
    randn('seed',sum(100*clock));
end

[nyear,nday]=size(seq);
sample0=seq(:,1:nday-1);
sample1=seq(:,2:nday);
nstate=max(max(seq));

for i=1:nstate
    for j=1:nstate
        NB(i,j)=length(find(sample0==i & sample1==j));
    end
end

NB=NB';
nb_tot=sum(NB);
NB=NB./(ones(nstate,1)*nb_tot);

NBC=cumsum(NB);
if nyear==1
    Y(1,1:nsim)=seq(1)*ones(1,nsim);
else
    Y(1,1:nsim)=randsample(seq(:,1)',nsim,'true');
end

for i=1:nsim; % loop on the number of simulated time series
    for j=2:nday;
        u=rand(1,nsim);
        tms=NBC(:,Y(j-1,i));
        a=sort([tms;u(i)]);
        b=find(a==u(i));
        Y(j,i)=b;
    end
end

write_ascii(out,Y);