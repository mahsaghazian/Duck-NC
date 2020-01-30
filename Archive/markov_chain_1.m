function [occ]=markov_chain_1(value,n,p01,p11,nsim,out);

% [occ]=markov_chain_1(value,n,p01,p11,nsim,out);
%
% Genegration of sequences of 0 and 1 that obeys to a Markov chain of order
% 1 (i.e. the state at t depends only on state at day at t-1). The 
% generation of sequences follows the method of Wilks (1998, Journal of 
% Hydrology, 178-191, eq. 3-5). 
%
% Inputs 
% 'value' : real number to initiate the random sequence (if > 0, the seeds
% is initiated to the value; otherwise, it is changed from the clock)
% 'n' : integer number = length of simulated time series 
% 'p01' and 'p11' : real number that are respectively the
% probability transition from state 0 to state 1 and the persistance of
% state 1 
% 'nsim' : integer number = number of time series you want
% to generate.
% 'out' : character string = name of the output file
%
% Outputs 
% 'occ' = matrice of n rows and nsim columns of 0 (= state 0) and 1 
% (= state 1) 
%
% Vincent MORON
% December 2004

if(value>0);
    rand('seed',value);
    randn('seed',value);
else
    rand('seed',sum(100*clock));
    randn('seed',sum(100*clock));
end

for i=1:nsim;
    a=rand(1); if a <= 0.5; occ(1,i)=0; else; occ(1,i)=1; end
    for j=2:n;
       if occ(j-1,i)==0; 
           Pcrit=p01; 
       else
           Pcrit=p11; 
       end;
       s=rand(1);
       if s <= Pcrit;
           occ(j,i)=1; 
       else 
           occ(j,i)=0; 
       end;
    end
end

if nargin==6
   write_ascii(out,occ);
end



