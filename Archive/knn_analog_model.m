function [TARGET_SIM,INDICEDAY]=knn_analog_model(ERA,MODEL,RAIN,prop,lseason,window,knn,value);

% [TARGET_SIM,INDICEDAY]=knn_analog_model(OBS,MODEL,TARGET,prop,lseason,window,knn,value);
%
% DRIVER for computing knn analog analysis between observed and simulated
% fields to estimate a target observed variable. 
% 
% Inputs 
% 'OBS' : a matrix of real number with rows describing time (typically
% daily values with seasons beneath each others) and columns describing
% space.
% 'MODEL' : a matrix of real number with rows describing time (typically
% daily values with seasons beneath each others) and columns describing
% space. The columns of 'MODEL' should refer to the same points as thos of
% 'OBS'. If there arre multiple runs, they should be placed beneath each
% others.
% 'TARGET' : a matrix of real number with rows describing time (typically
% daily values with seasons beneath each others) and columns describing
% space. The number of rows of 'TARGET' is the same as 'OBS'.
% 'prop' : scalar integer giving the proportion of variance explained by
% the EOF.
% 'lseason' : scalar integer giving the length of each season.
% 'window' : scalar integer giving the size of the moving temporal window
% 'knn' : scalar integer giving the number of simulation for each day.
% 'value' : real number to initiate the random sequence (if > 0, the seeds
% is initiated to the value; otherwise, it is changed from the clock)
%
% Outputs
% 'TARGET_SIM' : a matrix of real number giving the simulations  of the
% 'target' for each day. the first 'knn' rows are the simulations for the
% first day, then the next 'knn' rows are the simulations for the second
% day ...
% 'INDICEDAY' : a vector of integer giving the rank of each simulated day.
%
% Vincent MORON (dec 2004. version 1)
% modif. 5/1/2005 : knn neighbour are selected according to their
% relative position (the closest analog is selected the most frequently)
% modif. 12/1/2005 : verification of the steps : OK

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
R=R(PR);
PCS=PCS(:,PR)*diag(R);
PCSM=PCSM(:,PR)*diag(R);

PCPADD=padd(PC,lseason,window/2);
RAINPADD=padd(RAIN,lseason,window/2);
num=reshape([1:nyear*(lseason+window)],lseason+window,nyear);

TARGET_SIM=NaN*ones(nt2*knn,nv3);

% STEP 3 : searching for analogs in the sample of years
for i=1:nyear*nrun;
    disp(i);
    clear UNI pos target_model sample_mod posmod posmod_resample
    d=distance_euclid(PCS(:,:),PCSM(i,:),'sqeuclid');
    d2=(1./(d./max(d)));
    W2=cumsum(d2./sum(d2));
    UNI=rand(nyear,1);
    for j=1:nyear; % loop to resample years according to the ED with the target seasons
        [a,b]=sort([W2(:);UNI(j)]);
    	pos(j)=find(a==UNI(j));
        pos_resample(:,j)=num(:,pos(j));
    end
    target=PCM(((i-1)*lseason)+1:i*lseason,:);
    rainmod=NaN*ones(lseason*knn,nv3);
    for j=1:lseason;
        pos_resample2=pos_resample(j:j+window,:);
        library=PCPADD(pos_resample2,:);
        rain_library=RAINPADD(pos_resample2,:);
        [o,target_model(((j-1)*knn)+1:j*knn,:)]=find_analog(target(j,:),library,rain_library,knn,value); 
        indice_day2(((j-1)*knn)+1:j*knn)=pos_resample2(o);
    end
    TARGET_SIM(((i-1)*lseason*knn)+1:i*lseason*knn,1:nv3)=target_model;
    INDICEDAY(((i-1)*lseason*knn)+1:i*lseason*knn)=indice_day2;
end
