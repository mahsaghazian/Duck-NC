function [P,Q,HEP,HEQ,S,A,B]=bp(X,Y,choicestand,choice,PROP,latx,laty);

% [P,Q,HEP,HEQ,S,A,B]=bp(X,Y,choice1,choice2,latx,laty,PROP);
%
% CANONICAL CORRELATIONS ANALYSIS WITH EOF PREFILTERING
%
% Inputs ;
% 'X' and 'Y' are input matrices  of real numbers with
% columns usually describing time and rows describing space.
% 'choice1' : string character defining the standardisation of the matrices
% 'X' and 'Y' (='m' for normalization to zero mean and ='s' for
% standardization to zero mean and unit variance).
% 'choice2' : string character defining the weighting of the points with
% the latitudes (='weighting' if you want to weight each column of 'X' and
% 'Y' with the latitude)
% 'latx' : vector of real number giving the latitudes of 'X'
% 'laty' : vector of real number giving the latitudes of 'Y'
% 'PROP' : scalar indicating the proportion of variance you need to retain
% Bretherton et al (1992) considers a proportion of 0.7, but you need
% test the value which is the more accurate for your problem.
%
% Outputs ;
% 'P' contains the patterns (with columns describing order) of X
% which are also the homogenous covariance of X
% 'Q' contains the patterns (with columns describing order) of Y
% which are also the homogenous covariance of Y
% 'HEP' contains the heteregenous covariance of X
% 'HEQ' contains the heteregenous covariance of Y
% 'S' contains the canonical correlations (= singular value)
% 'A' contains the time expansion of U
% 'B' contains the time expansion of V
% The correlation between A and B contains the singular value on
% the main diagonal.
%
% see Intercomparison of methods for finding coupled patterns in climate
% data. Bretherton et al., 1992, J of Climate 5, 541-560.
%
% Vincent Moron (15/12/1998, checked 04/01/1999)
% moron@imga.bo.cnr.it

[Rx,Cx]=size(X);
[Ry,Cy]=size(Y);

if Rx ~= Ry
   error('X and Y must be of the same length')
end

X=stan(X,choicestand);
Y=stan(Y,choicestand);

Y=Y-(ones(Ry,1)*mean(Y));
if strcmp(choice,'weighting');
    X=ponder(X,latx);
    Y=ponder(Y,laty);
end

% decomposition of both matrix into EOFs
[Vx,ux]=eig(cov(X)); [ii,jj]=sort(diag(ux)); Vx=fliplr(Vx(:,jj)); ux=flipud(ii);
[Vy,uy]=eig(cov(Y)); [ii,jj]=sort(diag(uy)); Vy=fliplr(Vy(:,jj)); uy=flipud(ii); 
indx=(find(cumsum(ux) > PROP*sum(ux)));
indy=(find(cumsum(uy) > PROP*sum(uy)));
px=(sum(ux(1:indx(1)))/sum(ux))*100;
py=(sum(uy(1:indy(1)))/sum(uy))*100;
disp([num2str(indx(1)),' EOFs are retained for X explaining ',num2str(px),' of the variance'])
disp([num2str(indy(1)),' EOFs are retained for Y explaining ',num2str(py),' of the variance'])
ux=diag(ux(1:indx(1)));
uy=diag(uy(1:indy(1))); 
Vx=Vx(:,1:length(ux));
Vy=Vy(:,1:length(uy));
Tx=X*Vx;
Ty=Y*Vy;
Tx=Tx-(ones(Rx,1)*mean(Tx));
Ty=Ty-(ones(Ry,1)*mean(Ty));
Tx=Tx/sqrt(ux); % normalization of PCs of X
Ty=Ty/sqrt(uy); % normalization of PCs of Y

% computation of the cross-covariance matrices
Cxy=(Tx'*Ty)./(Rx-1);
[U,S,V]=svd(Cxy); % S contains canonical correlations on the main diagonal

% patterns = homogenous covariance
P=Vx*(ux^0.5)*U;
Q=Vy*(uy^0.5)*V;

% time expansion coefficients
A=Tx*U;% time expansion for the left field
B=Ty*V;% time expansion for the right field

% heteregenous covariance
mini=min([size(P),size(Q)]);
S=diag(S);
S=S(1:mini);
HEP=P(:,1:mini)*diag(S);
HEQ=Q(:,1:mini)*diag(S);
