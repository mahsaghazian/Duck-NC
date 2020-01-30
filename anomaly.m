function [Y,CYCM,CYCS]=anomaly(X,lseason,nrun,cpu,choice)

% [Y,CYCM,CYCS]=anomaly(X,lseason,nrun,cpu,choice)
%
% This function computes anomaly of daily data (or instataneous data with a
% daily timestep. 
% 
% Input :
% 'X': matrix containing data with column describing space and rows 
% describing time (the dimensions are ordered as day,year and run with the 
% first one varying the fastest). 
% 'lseason' : integer giving the length of the season,
% 'nrun' : integer is the number of runs (=1 in case of observations or 
% ERA40). If there are mutliple runs, each one is processed
% independently to the others. 
% 'cpu' : real scalar giving the cut-off frequency to filter the seasonal 
% cycle (th daily mean is firstly computed on all year available
% and then filtered). 
% 'choice' : character string that could be 'm' or 's'; If it is 'm', the
% anomaly are computed respetively to the filtered seasonal cycle and if it
% is 's', the anomaly are divided by a filtered seasonal cycle of the
% standard deviations (computed as the seasonal cycle of the mean).
% 
% Output :
% 'Y' : matrix containing the anomalies
% "CYCM' : matrix containing the mean seasonal cycle of the mean values
% 'CYCS' : matrix containing the mean seasonal cycle of the std values
% (empty if choice='m');
%
% created by; Vincent Moron
% August 2005

[Nt,Nv]=size(X);
Nyear=Nt/(lseason*nrun);
Nt2=lseason*Nyear;

for i=1:nrun;
    sample=X(((i-1)*Nt2)+1:i*Nt2,:);
    for j=1:lseason;
        m(j,:)=mean(sample(j:lseason:Nt2,:));
        s(j,:)=std(sample(j:lseason:Nt2,:));
    end
    bfm=filtrage(m,9,cpu,'low');
    bfs=filtrage(s,9,cpu,'low');
    bfm=copy(bfm,Nyear);
    bfs=copy(bfs,Nyear);
    if strcmp(choice,'m');
        Y(((i-1)*Nt2)+1:i*Nt2,:)=sample-bfm;
        CYCM(((i-1)*Nt2)+1:i*Nt2,:)=bfm;
        CYCS=[];
    elseif strcmp(choice,'s');
        Y(((i-1)*Nt2)+1:i*Nt2,:)=s(sample-bfm)./bfs;  
        CYCM(((i-1)*Nt2)+1:i*Nt2,:)=bfm;
        CYCS(((i-1)*Nt2)+1:i*Nt2,:)=bfs;
    end
end