function [CLUS,mean_cluster,CLUSMOD,RAINMOD,ID]=weather_type_model(ERA,MODEL,RAIN,C,CM,CMOD,missing,prop,nc,lseason,window,knn,value);

% [clus,mean_cluster,clusmod,TARGET_SIM,ID]=weather_type_model(OBS,MODEL,TARGE
% T,RAIN,prop,nc,lseason,window,knn,value,C,CM,CMOD);
%
% DRIVER for computing weather type classification between observed and 
% simulated fields to estimate a target observed variable. 
% 
% Inputs 
% 'OBS' : a matrix of real number with rows describing time (typically
% daily values with seasons beneath each others) and columns describing
% space.
% 'MODEL' : a matrix of real number with rows describing time (typically
% daily values with seasons beneath each others) and columns describing
% space. The columns of 'MODEL' should refer to the same points as those of
% 'OBS'. If there arre multiple runs, they should be placed beneath each
% others.
% 'TARGET' : a matrix of real number with rows describing time (typically
% daily values with seasons beneath each others) and columns describing
% space. The number of rows of 'TARGET' is the same as 'OBS'.
% 'prop' : scalar integer giving the proportion of variance explained by
% the EOF.
% 'nc' : number of clusters
% 'lseason' : scalar integer giving the length of each season.
% 'window' : scalar integer giving the size of the moving temporal window
% 'knn' : scalar integer giving the number of simulation for each day.
% 'value' : real number to initiate the random sequence (if > 0, the seeds
% is initiated to the value; otherwise, it is changed from the clock)
% 'C' : vector of integer giving the identifier of clusters in 'OBS'
% 'CM' : vector of integer giving the identifier of cluster in 'MODEL'
% 'CMOD' : matrix of real number (= centroid of clusters in the EOF-space
% of 'OBS'). If 'C', 'CM' and 'CMOD' are provided as empty scalar (=[])
% they are computed using 'OBS' and parameters defined above). 
%
% Outputs
% 'clus' : vector of integer giving the identifier of clusters in 'OBS'
% 'mean_cluster' : matrix of real number (= centroid of clusters in the EOF-space
% 'clusmod' : vector of integer giving the identifier of cluster in 'MODEL'
% 'TARGET_SIM' : a matrix of real number giving the simulations  of the
% 'target' for each day. the first 'knn' rows are the simulations for the
% first day, then the next 'knn' rows are the simulations for the second
% day ...
% 'ID' : a vector of integer giving the rank of each simulated day.
%
% Vincent MORON (dec 2004. version 1)

if(value>0);
    rand('seed',value);
    randn('seed',value);
else
    rand('seed',sum(100*clock));
    randn('seed',sum(100*clock));
end

[nt1,nv1]=size(ERA);
[nt2,nv2]=size(MODEL);
[nt3,nv3]=size(RAIN);

nyear=nt1/lseason;
nrun=nt2/(lseason*nyear);
RAINMOD=NaN*ones(nt2*knn,nv3);

% standardisation of matrices and computation of seasonal mean
ERAS=stan(ERA,'s');
MODELS=NaN*ones(size(MODEL));
for i=1:nrun;
    MODELS(((i-1)*nt1)+1:i*nt1,:)=stan(MODEL(((i-1)*nt1)+1:i*nt1,:),'s');
end

ERAM=reshape(mean(reshape(ERA,lseason,nyear*nv1)),nyear,nv1);
MODELM=reshape(mean(reshape(MODEL,lseason,nyear*nrun*nv2)),nyear*nrun,nv2);
ERAM=stan(ERAM,'s');
for i=1:nrun;
    MODELM(((i-1)*nyear)+1:i*nyear,:)=stan(MODELM(((i-1)*nyear)+1:i*nyear,:),'s');
end
    
% pre-filtering with EOF on daily value    
[U,S,V]=svd(ERAS,0);
SC=diag(S).^2; 
SC=SC./sum(SC);
a=find(cumsum(SC) > prop); A=a(1);
disp([num2str(A),' daily PCs represent ',num2str(prop*100),'% of the variance']);
PC=U(:,1:A)*S(1:A,1:A);
V=V(:,1:A);
std_pc=std(PC);
PCM=MODELS*V; % simulated data are projected onto the EOF of ERA

% the standard deviation of PCM is scaled
for i=1:nrun;
    ratio=std_pc./std(PCM(((i-1)*lseason*nyear)+1:i*lseason*nyear,:));
    ratio=diag(ratio);
    PCM(((i-1)*lseason*nyear)+1:i*lseason*nyear,:)=PCM(((i-1)*lseason*nyear)+1:i*lseason*nyear,:)*ratio;
end

% pre-filtering with EOF on seasonal values
[UM,SM,VM]=svd(ERAM,0);
SCM=diag(SM).^2; 
SCM=SCM./sum(SCM);
am=find(cumsum(SCM) > prop); AM=am(1);
disp([num2str(AM),' seasonal PCs represent ',num2str(prop*100),'% of the variance']);
PCS=UM(:,1:AM)*SM(1:AM,1:AM);
VS=VM(:,1:AM);
std_pcs=std(PCS);
PCSM=MODELM*VS;

for i=1:nrun;
    ratios=std_pcs./std(PCSM(((i-1)*nyear)+1:i*nyear,:));
    ratios=diag(ratios);
    PCSM(((i-1)*nyear)+1:i*nyear,:)=PCSM(((i-1)*nyear)+1:i*nyear,:)*ratios;
end

% skill of seasonal PC
for i=1:AM;
    if nrun > 1
       R(i)=(1/(nyear-1))*stan(PCS(:,i),'s')'*stan(mean(reshape(PCSM(:,i),nyear,nrun)'),'s')';
    else
       R(i)=(1/(nyear-1))*stan(PCS(:,i),'s')'*stan(PCSM(:,i),'s');
    end
end

PR=find(R > 0);
disp([num2str(length(PR)),' seasonal PCs, having skill > 0, are retained']);
R=R(PR)
PCS=PCS(:,PR)*diag(R);
PCSM=PCSM(:,PR)*diag(R);

% k-means classification with 1000 iterations
if C(1)==missing;
    [clus,mean_cluster]=kmeans(PC,nc,'Replicates',10,'Maxiter',10000);
else
    clus=C;
    mean_cluster=CM;
end

% computation of the centroid of each cluster on the data
disp('K-means classification done ...');

if CMOD(1)==missing;
    for i=1:lseason*nyear*nrun;
        d=distance_euclid(mean_cluster,PCM(i,:),'sqeuclid');
        clusmod(i)=find(d==min(d));
    end
else
    clusmod=CMOD;
end
disp('GCM day are classified ...')

CLUS=clus;
CLUSMOD=clusmod;

clus=padd(clus,lseason,window/2);
RAINPADD=padd(RAIN,lseason,window/2);
na=0;
nb=0;

clus=reshape(clus,lseason+window,nyear);
clusmod=clusmod(:);
num=reshape([1:nyear*(lseason+window)],lseason+window,nyear);

for i=1:nyear*nrun;
    disp(i);
    clear UNI pos clus_resample pos_rain sample_mod posmod posmod_resample
    d=distance_euclid(PCS(:,:),PCSM(i,:),'sqeuclid');
    d2=(1./(d./max(d)));
    W2=cumsum(d2./sum(d2));
    UNI=rand(100,1);
    for j=1:100; % loop to resample years according to the ED with the target seasons
        [a,b]=sort([W2(:);UNI(j)]);
    	pos(j)=find(a==UNI(j));
        clus_resample(:,j)=clus(:,pos(j));
        pos_rain(:,j)=num(:,pos(j));
    end
    sample_mod=clusmod(((i-1)*lseason)+1:i*lseason);
    rainmod=NaN*ones(lseason*knn,nv3);
    clear j
    for j=1:lseason;
        clus_resample2=clus_resample(j:j+window,1:nyear);
        clus_resample2=clus_resample2(:);
        pos_rain2=pos_rain(j:j+window,:);
        pos2=find(clus_resample2==sample_mod(j));
        if isempty(pos2);% relaxation to the whole seasson if pos2 is empty
            na=na+1;
            clus_resample2=clus_resample(:);
            pos_rain2=pos_rain(:);      
            pos2=find(clus_resample2==sample_mod(j));
            if isempty(pos2); % relaxation to whole period if pos2 is still empty
                nb=nb+1;
                pos2=find(clus==sample_mod(j));
                pos2=randsample(pos2,knn,'true');
                pos_rain2=num;
            else
                pos2=randsample(pos2,knn,'true');
            end
        else
            pos2=randsample(pos2,knn,'true');
        end
        rainmod(((j-1)*knn)+1:j*knn,:)=RAINPADD(pos_rain2(pos2),:);
        id2(((j-1)*knn)+1:j*knn)=(pos_rain2(pos2));
    end
    RAINMOD(((i-1)*lseason*knn)+1:i*lseason*knn,:)=rainmod;
    ID(((i-1)*lseason*knn)+1:i*lseason*knn)=id2;
end                 
disp(['the search is relaxed to the whole season in ',num2str(na),' cases']);
disp(['the search is relaxed to the whole sample in ',num2str(nb),' cases']);

