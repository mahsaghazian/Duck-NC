function [A_verif,B_verif,ADJ,ADJM]=mos_cca(X,Y,stand,choice,latx,laty,DIMS,PROP);

% [A,B,ADJ,ADJM]=mos_cca(X,Y,stand,choice,latx,laty,DIMS,PROP);
%
% Cross validated reconstruction of a predictand field 'X' from a predictor
% filed 'Y'. 'X' could be observations and 'Y' could be run of ensembles
% (they should be placed beneath each others in that case) but this
% function works also with two observed fields. The reconstruction is based
% on a CCA between both fields.
%
% Input
% 'X' : matrix of real numbers describing the left field (NOT STANDARDIZED)
% 'Y' : matrix of real numbers describing the right field (NOT 
% STANDARDIZED). If 'Y' contains multiple runs, they need to be beneath
% each others.
% 'stand' : string character defining the standardisation of the matrices
% 'X' and 'Y' (='m' for normalization to zero mean and ='s' for
% standardization to zero mean and unit variance). In case of multiple runs
% they are standardized independently.
% 'choice' : string character defining the weighting of the points with
% the latitudes (='weighting' if you want to weight each column of 'X' and
% 'Y' with the latitude)
% 'latx' : vector of real number giving the latitudes of 'X'
% 'laty' : vector of real number giving the latitudes of 'Y'
% 'DIMS' : integer giving the dimension in unit time of the verification
% periods.
% 'PROP' : scalar indicating the proportion of variance you need to retain
% in the EOF prefiltering.
%
% Outout
% 'A' : matrix of real number giving the cross-validated time scores for
% 'X'
% 'B' : matrix of real number giving the cross-validated time scores for
% 'Y'
% 'ADJ' : matrix of real number giving the MOS reconstructions of 'X' using
% 'Y' as a predictor. The first nc columns (nc = number of variables of
% 'X') is the reconstructions using the leading CCA mode, then the next nc
% columns are the reconstructions using the second CCA mode and so on. The
% sum of the reconstructions gives the MOS for the all CCA modes, but in
% general, it is better to use the CCA modes that have a positive (or
% significant positive) correlation between 'A' and 'B'.
% 'ADJM' : same as 'A' except that for the mean of the ensemble. In case of
% two observed fields, ADJM is ADJ.
%
% Vincent Moron
% March 1997

[nyear,varx]=size(X);
[Ry,vary]=size(Y);
nrun=Ry/nyear;
disp(['there are ',num2str(nrun),' runs...']);

XX=stan(X,stand); 
XX=copy(XX,nrun);
YY=stan(reshape(Y,nyear,nrun*vary),stand); 
YY=reshape(YY,nyear*nrun,vary);

YY(find(isinf(YY) | isnan(YY)))=zeros(size(find(isinf(YY) | isnan(YY))));
XX(find(isinf(XX) | isnan(XX)))=zeros(size(find(isinf(XX) | isnan(XX))));

[u,v,hep,heq,s,a,b]=bp(XX,YY,stand,choice,PROP,latx,laty);
[nx,indx]=size(u);
[ny,indy]=size(v);
ls=min([indx indy]);
disp('CCA whitout cross-validation done ...');
disp(['there are ',num2str(ls),' reconstructions']);

% definition of the cross-validated period
[CV]=cv(nyear,DIMS);
[nl,nper]=size(CV);
CV=copy(CV,nrun);

disp(['there are ',num2str(nper),' periods of cross-validation']);

% All the cross-validated periods 
for i=1:nper
   	ind_training1=find(CV(1:nyear,i)==0); l_t1=length(ind_training1);
    ind_training2=find(CV(:,i)==0); l_t2=length(ind_training2);
    ind_verif1=find(CV(1:nyear,i)==1); l_v1=length(ind_verif1);
    ind_verif2=find(CV(:,i)==1); l_v2=length(ind_verif2);
    X_training=X(ind_training1,:);
    Y_training=Y(ind_training2,:);
    MY_training=reshape(copy(mean(reshape(Y_training,l_t1,vary*nrun)),l_v1),l_v1*nrun,vary);
    SY_training=reshape(copy(std(reshape(Y_training,l_t1,vary*nrun)),l_v1),l_v1*nrun,vary);
    X_verif=X(ind_verif1,:);
    Y_verif=Y(ind_verif2,:);
    % standardization of the matrices ; each run is processed independently
    XS_training=stan(X_training,stand); 
    XS_training=copy(XS_training,nrun);
    if strcmp(choice,'weighting'); [XS_training,wx]=ponder(XS_training,latx); end 
    YS_training=reshape(stan(reshape(Y_training,l_t1,nrun*vary),stand),l_t1*nrun,vary);
    if strcmp(stand,'m');
          XS_verif=X_verif-(ones(l_v1,1)*mean(X_training));
    else
          XS_verif=(X_verif-(ones(l_v1,1)*mean(X_training)))./(ones(l_v1,1)*std(X_training));
    end    
    XS_verif=copy(XS_verif,nrun);
    if strcmp(choice,'weighting'); [XS_verif,wx]=ponder(XS_verif,latx); end 
    if strcmp(stand,'m');
          YS_verif=Y_verif-MY_training;
    else
          YS_verif=(Y_verif-MY_training)./(SY_training);
    end   
    YS_training(find(isinf(YS_training) | isnan(YS_training)))=zeros(size(find(isinf(YS_training) | isnan(YS_training))));   
    XS_training(find(isinf(XS_training) | isnan(XS_training)))=zeros(size(find(isinf(XS_training) | isnan(XS_training))));
    YS_verif(find(isinf(YS_verif) | isnan(YS_verif)))=zeros(size(find(isinf(YS_verif) | isnan(YS_verif))));   
    XS_verif(find(isinf(XS_verif) | isnan(XS_verif)))=zeros(size(find(isinf(XS_verif) | isnan(XS_verif))));
    if strcmp(choice,'weighting'); [YS_verif,wx]=ponder(YS_verif,laty); end 
    if strcmp(choice,'weighting'); [YS_training,wx]=ponder(YS_training,laty); end 
    % CCA on each cross-validated period
    [Vx,ux]=eig(cov(XS_training)); [ii,jj]=sort(diag(ux)); Vx=fliplr(Vx(:,jj)); ux=flipud(ii);
    [Vy,uy]=eig(cov(YS_training)); [ii,jj]=sort(diag(uy)); Vy=fliplr(Vy(:,jj)); uy=flipud(ii); 
    ux=diag(ux(1:indx));
    uy=diag(uy(1:indy));
    Vx=Vx(:,1:length(ux));
    Vy=Vy(:,1:length(uy));
    Tx=XS_training*Vx;
    Ty=YS_training*Vy;
    Txx=XS_verif*Vx;
	Tyy=YS_verif*Vy;
	Tx=Tx-(ones(l_t2,1)*mean(Tx));
	Ty=Ty-(ones(l_t2,1)*mean(Ty));
	Tx=Tx/sqrt(ux); % normalization of PCs of X
	Ty=Ty/sqrt(uy); % normalization of PCs of Y
	Txx=Txx-(ones(l_v2,1)*mean(Tx));
    Tyy=Tyy-(ones(l_v2,1)*mean(Ty));
    Txx=Txx/sqrt(ux);
    Tyy=Tyy/sqrt(uy);	
% computation of the cross-covariance matrices
	Cxy=(Tx'*Ty)./(l_t2-1);
	[U,S,V]=svd(Cxy); % S contains canonical correlations on the main diagonal
    UU=Vx*(ux^0.5)*U; UU=UU(:,1:ls);
	VV=Vy*(uy^0.5)*V; VV=VV(:,1:ls);   
    P(:,((i-1)*ls)+1:i*ls)=UU;
    Q(:,((i-1)*ls)+1:i*ls)=VV;
    U_training(:,((i-1)*ls)+1:i*ls)=U(:,1:ls);
    V_training(:,((i-1)*ls)+1:i*ls)=V(:,1:ls);
    A_training=Tx*U_training(:,((i-1)*ls)+1:i*ls);
    B_training=Ty*V_training(:,((i-1)*ls)+1:i*ls);
    Txx_verif(ind_verif2,:)=Txx(:,1:indx);
    Tyy_verif(ind_verif2,:)=Tyy(:,1:indy);     
    U_eof=U_training(:,((i-1)*ls)+1:i*ls);
    V_eof=V_training(:,((i-1)*ls)+1:i*ls);
    A_verif(ind_verif2,:)=Txx_verif(ind_verif2,:)*U_eof;
    B_verif(ind_verif2,:)=Tyy_verif(ind_verif2,:)*V_eof;
    S_training(:,i)=diag(S);
    P_training=P(:,((i-1)*ls)+1:i*ls);
    for j=1:ls,
        ADJ(ind_verif2,((j-1)*varx)+1:j*varx)=B_verif(ind_verif2,j)*S(j,j)*P_training(:,j)';
    end
end

% re-ordonnement 
[ra,ca]=size(A_verif);
R=(1/(ra-1))*stan(A_verif,'s')'*stan(B_verif,'s');
[R,order]=sort(diag(R));
order=flipud(order);
A_verif=A_verif(:,order);
B_verif=B_verif(:,order);
[radj,cadj]=size(ADJ);
ADJ=reshape(ADJ,radj*varx,cadj/varx);
ADJ=ADJ(:,order);
ADJ=reshape(ADJ,radj,cadj);

if nrun > 1
    ADJM=reshape(mean(reshape(reshape(ADJ,nyear,nrun*varx*ls)',nrun,nyear*varx*ls)),varx*ls,nyear);
    ADJM=ADJM';
else
    ADJM=ADJ;
end
  


