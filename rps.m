function [RPS,RPSS]=rps(obs,sim,nc);

% [RPS,RPSS]=rps(obs,sim,nc);
%
% This function compute the rank probability skill score. This score is a
% measure of skill between an observed vector and an ensemble of
% simulations. The score is based on the cumulative sum of categorical
% forecasts.
% 
% Input
% 'obs' : a vector of real number describing the temporal evolution of any
% variable
% 'sim' : a matrix of real number where the columns describes the same time
% unit as in 'obs' and columns give the different estimate for each time
% increment. Note that 'obs' and 'sim' should be scaled in the same units.
% 'nc' : a vector giving the number of observations from the lowest to the 
% highest value, in each class (sum(nc) is the number of rows of 'obs' and
% 'sim');
%
% Output
% 'RPS' : a vector of real number giving rank probability score (= 0 if all
% the simulations are in the same class as observations)
% 'RPSS' : a vector of real number giving rank probability skill score for
% each time unit (= 100% if RPS = 0 ; = 0% if the skill is the same as a 
% climatological forecast and < 0% ifthe forecast is worst than a 
% climatological forecast).
% 
% Vincent MORON
% Fev. 2004

obs=obs(:);
[nt,nr]=size(sim);
sim=stan(sim,'s');
obs=stan(obs,'s');
[a,b]=sort(obs);
ncs=cumsum(nc);

% �tape 1 : classement des observations et calcul des seuils
for l=1:length(nc)-1;
    thres(l)=a(ncs(l)+1); % 'thres' contient les seuils de la classification
end
thres=[-Inf thres +Inf];
clobs=NaN*ones(size(obs));
for i=1:length(nc);
    aa=find(obs >= thres(i) & obs < thres(i+1));
    clobs(aa)=i*ones(size(aa));
end

% �tape 2 : cr�ation de la matrice du classement des observations
classo=zeros(nt,length(nc));
for j=1:nt;
    classo(j,clobs(j))=[1]; % cette matrice contient un 1 et des 0 par ligne
end

% �tape 3 : classement des simulations PDF
for j=1:nt;
    simul=sim(j,:);
    for i=1:length(nc);
         aa=find(simul >= thres(i) & simul < thres(i+1));
         cl(i)=length(aa);
    end
    clsim(j,:)=cl; clear cl aa simul
end

clsim=clsim./nr;

% �tape 4 : calcul des PDFs cumul�es et du RPS
classo=cumsum(classo');
clsim=cumsum(clsim');
RPS=sum((clsim-classo).^2);

% �tape 5 : calcul de la climatologie et de la RPS climatologique
climo=nc./nt;
climo=cumsum(climo(:));
climo=climo*ones(1,nt);
RPSclim=sum((climo-classo).^2);
RPSS=100*(1-(RPS./RPSclim));

