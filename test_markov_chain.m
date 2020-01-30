function [prob,plike,punlike]=test_markov_chain(seq,value,nsim);

% [prob,plike,punlike]=test_markov_chain(seq,value,nsim)
%
% This prg tests the probability transition between states (identified
% through a scalar). It has been developed for daily time scales and the
% amtrix "seq" contains the states with row describing years and column
% descring days. nsim is the number of simulations needed. The output prob
% gives the observed probability from state i (in row) to state j (in
% column) and p is the monte carlo significance. The significance is
% determined using a permuted matrix with the same frequency for each state
% as observed (Vautard et al., J. Atmos. Sci., 1990). 
%
% Inputs
% 'seq' : input vector/matrix to be tested. The rows describe the years
% and the colums describes the days and the matrix contains integers.
% 'value' : real number to initiate the random sequence (if > 0, the seeds
% is initiated to the value; otherwise, it is changed from the clock)
% 'nsim' : integer number giving the number of simulations
% 
% Outputs
% 'prob' : matrix of real number giving the probability transition
% 'plike': matrix of the probability for which the transition from state i 
% to state j in the  random ensemble is below the observed value
% 'punlike' : matrix of the  probability for which the transition from 
% state i to state j in the random ensemble is above the observed value. 
% Vincent Moron
% feb. 2005

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
        nb(i,j)=length(find(sample0==i & sample1==j));
    end
end

nb_tot=sum(nb);
prob=nb./(ones(nstate,1)*nb_tot);
nomb=round((nyear*nday)*(nb_tot/sum(nb_tot)));
cnomb=[0,cumsum(nomb)];

seqsim=NaN*ones(1,nyear*nday);
for i=1:nstate
    ll=(cnomb(i+1)-cnomb(i));
    seqsim(cnomb(i)+1:cnomb(i+1))=i*ones(1,ll);
end

for k=1:nsim;
    r=seqsim(randperm(nday*nyear));
    r=r(1:nday*nyear);
    r=reshape(r,nyear,nday);
    sample0=r(:,1:nday-1);
    sample1=r(:,2:nday);
    for i=1:nstate
        for j=1:nstate
            nb_sim(i,j)=length(find(sample0==i & sample1==j));
        end
    end
    NB(:,k)=nb_sim(:);
end

nbobs=nb(:);
for i=1:length(nbobs);
    plike(i)=length(find(NB(i,:) >= nbobs(i)))/nsim;
    punlike(i)=length(find(NB(i,:) <= nbobs(i)))/nsim;
end

plike=reshape(plike,nstate,nstate);
punlike=reshape(punlike,nstate,nstate);
