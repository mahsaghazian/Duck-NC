function [hr,far,rocscore]=roc_prob(obs,forecast,ncl,proba);

% [HR,FAR,ROCSS]=rocscore(obs,forecast,ncl,prob) 
%
% This function computes ROC scores (hit rate, false-alarm rate and roc
% score) for probabilistic forecasts. the observed data are firtsly 
% classified into 'ncl' equiprobable classes. Then runs are distributed 
% into these classes and the percentage of runs falling into each classe is 
% computed for each time step. Hit Rate (HR), and False-Alarm (FAR) rates 
% are computed for each probability of warning. the ROC score is simply = 
% 2 * (A-0.5) where A is the surface below the ROC curve. Ths score = 1 
% for a perfect forecast, 0 for no skill and -1 for a perfectly bad 
% forecast.
%
% Input 
% 'obs' : a vector of real number giving the observed value
% 'forecast' : a matrice of real number giving the estimated value of the
% observations with rows describing time (= length of vector 'obs') and
% columns describing the different estimates (as different runs for
% example)
% 'ncl' : number of classes (= 3 for having three classes
% for example. The classes are equiprobable and their threshold are 
% computed from percentiles)
% 'prob' : probability on which hit and false-alarm rates are
% computed (in percentages as [0:10:100]) 
%
% Output
% 'HR' : vector of real number giving the hit-rate for each 'prob'
% 'FAR' : vector of real number giving the False-alarm rate for each
% 'prob'
% 'rocscore' : integer scalar giving the roc score
%
% Vincent MORON
% May 2005

warning off MATLAB:divideByZero
obs=obs(:);
[nyear,nrun]=size(forecast);
threshold=[1/ncl:1/ncl:1-(1/ncl)];
threshold=threshold*100;
lprob=length(proba);
[classe_obs,nobs,fobs]=classe(obs,prctile(obs,threshold));
[classe_for,nfor,ffor]=classe(forecast,prctile(obs,threshold));
classe_for=reshape(classe_for,nyear,nrun);

% computation of the probabilty of runs into each classe defined
% from observations. 

for i=1:ncl,
    for j=1:nyear,
        prob_sim(j,i)=(length(find(classe_for(j,:)==i))/nrun)*100;
    end
end

proba=flipud(proba(:));

for i=1:ncl,
    classe1=zeros(size(obs));
    classe1(find(classe_obs==i))=ones(size(find(classe_obs==i)));
    for j=1:lprob,
       classe2=zeros(size(obs));
       classe2(find(prob_sim(:,i) >= proba(j)))=ones(size(find(prob_sim(:,i) >= proba(j))));
       [hr(j,i),far(j,i),ll,tab]=roc(classe1,classe2);    
    end
    rocscore(i)=polyarea([0;far(:,i);1;1;0],[0;hr(:,i);1;0;0]);
end

rocscore=2.*(rocscore-0.5);


        
