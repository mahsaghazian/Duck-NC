function [GCM3,PGCM,SCALING]=local_scaling(gcm,obs,seas,threshold);

% [GCM2,Pgcm,scaling]=local_scaling(gcm,obs,seas,threshold);
%
% This function helps to scale/calibrate a GCM time series so that the
% frequency of wet days and mean intensity of rain during wet days matches
% the observed one. 
%
% Inputs
% 'gcm' : real matrix containing the simulated time series (if there are
% multiple runs, they need to be copied beneath each others)
% 'obs' : real matrix containing the observed time series
% 'seas' : scalar vector describing the months (for example, the July-
% September season should be described as 31 ones, then 31 two and then 30
% three
% 'threshold' : threshold for defining wet days
%
% 'Outputs'
% 'GCM2' : real vector containing the calibrated simulated time series
% 'Pgcm' : vector giving the threshold for each month
% 'scaling' : vector giving the scaling factor for each month
%
% Vincent Moron
% October 2005
% revised in March 2006 (for considering monthly-decadal scaling and a
% matrix instead of a vector)

[nt,nc]=size(obs);
[nt2,nc2]=size(gcm);
ls=length(seas);
nyear=nt/ls;
nrun=nt2/nt;
lm=max(seas);
S=copy(seas,nyear);
GCM2=NaN*ones(size(obs));

for i=1:nrun;
    disp(i);
    gcm2=gcm(((i-1)*nt)+1:i*nt,:);
    for j=1:nc;
        for k=1:lm;
            clear O G gcmr mobs mgcm GCM
            O=obs(find(S==k),j); % scaling on each months
            G=gcm2(find(S==k),j);
            l=length(find(O > threshold));
            gcmr=sort(G);
            gcmr=flipud(gcmr);
            Pgcm(k,j)=gcmr(l+1);
            % local scaling factor
            mobs=O(find(O > threshold));
            mobs=mean(mobs);
            mgcm=G(find(G > Pgcm(k,j)));
            mgcm=mean(mgcm);
            scaling(k,j)=(mobs-threshold)/(mgcm-Pgcm(k,j));
            GCM=(scaling(k,j)*(G-Pgcm(k,j)))+threshold;
            GCM(find(GCM < threshold))=zeros(size(find(GCM < threshold)));
            GCM2(find(S==k),j)=GCM(:);
        end
    end
    GCM3(((i-1)*nt)+1:i*nt,:)=GCM2;
    SCALING(((i-1)*lm)+1:i*lm,:)=scaling;
    PGCM(((i-1)*lm)+1:i*lm,:)=Pgcm;
    clear GCM2 scaling Pgcm;
end
