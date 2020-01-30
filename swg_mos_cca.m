function [RAINMOD]=swg_mos_cca(rain,mat,lseason,thres,latx,laty,cw,cdib,cs,DIMS,PROP,knn);

% [RAINMOD]=swg_mos_cca(rain,mat,lseason,thres,latrain,latmat,cw,cdib,cs,DIMS,PROP,knn);
%
% driver for generating daily sequences using a stochastic weather
% generator (one-order markov chain and 3 PDFs (gamma, Weibull and mixed
% exponential for the PDF of the rainy days). The statistic model use
% seasonal or monthly mean to predict the needed parameters (i.e.
% probability from dry to wet, probability from wet to wet, and the
% parameters of the PDF of rainy days) and then a SWG to generate daily
% sequences using these estimated values. The estimation is performed using
% a cross-validated CCA on monthly/seasonal values of 'rain' and 'mat'. The
% cross-validated CCA are independently performed on each parameter. In
% case of multiple runs in 'mat', the CCA is performed on the whole
% ensemble but estimates are computed using the mean ensemble. 
%
% Input
% 'rain' : matrix of real number giving the daily rainfall on a network of
% station (row = time and column = stations). Seasons need to be beneath
% each others.
% 'mat' : a matrix of real number giving the predictors of the rainfall
% parameters. It should be seasonal/monthly mean and in case of multiple
% runs, they should be placed beneath each others.
% 'lseason' : scalar integer giving the length of the season
% 'thres' : scalar integer giving the threshold of the definition of rainy
% days in 'rain' (i.e. a daily amount > threshold is wet)
% 'latrain' : vector of real number giving the latitudes of the stations
% in 'rain' if the CCA is performed on weighted values
% 'latmat' : vector of real number giving the latitudes of the predictors
% if the CCA is performed on weighted values.
% 'cw' : character string to define if CCA is performed on weighted values
% according to the latitudes of the points (='weighting' for a weighted CCA
% and 'no if not)
% 'cdib' : character string to define the PDF used to model amount during
% rainy days (='gamma', 'weibull' or 'mixed_exp')
% 'cs' : character string to define standardisation of data (='m' for 
% normalization to zero mean or 's' for standardization to zero mean and  
% unit variance.
% 'DIMS' : scalar integer to define the length of verification period of
% the cross-validated CCA
% 'PROP' : real giving the proportion of variance used to pre-filter the
% fields before the CCA with an EOF analysis.
% 'knn' : scalar giving the number of simulation for each season.
%
% Output
% 'RAINMOD' : matrix of real number giving the simulated rainfall (the
% 'knn' days are placed beneath each others and the columns have the same
% sense as in 'rain'.
%
% Vincent Moron
% Nov 2005

[nr,nc]=size(rain);
nyear=nr/lseason;
[nr2,nc2]=size(mat);
nrun=nr2./nyear

for i=1:nyear;
    [p00(i,:),p01(i,:),p10(i,:),p11(i,:)]=prob_wet_dry(rain(((i-1)*lseason)+1:i*lseason,:),thres);
end

for i=1:nyear;
    for j=1:nc;
        sample=rain(((i-1)*lseason)+1:i*lseason,j);
        sample=sample(find(sample > thres));
        if strcmp(cdib,'gamma');
            [g]=gamfit(sample);
            ga(i,j)=g(1);
            gb(i,j)=g(2);
            gc(i,j)=0;
        elseif strcmp(cdib,'weibull');
            [g]=wblfit(sample);
            ga(i,j)=g(1);
            gb(i,j)=g(2);
            gc(i,j)=0;
        elseif strcmp(cdib,'mixed_exp');
            [g]=mixed_expfit(sample,1);
            ga(i,j)=g(1);
            gb(i,j)=g(2);
            gc(i,j)=g(3);
        end
    end
end

p01z=log_odd(p01);
p11z=log_odd(p11);
for i=1:nc;
    a=p11z(:,i);
    b=find(isinf(a));
    c=find(~isinf(a));
    if ~isempty(b);
        p11z(b,i)=min(p11z(c,i));
    end
end
for i=1:nc;
    a=p01z(:,i);
    b=find(isinf(a));
    c=find(~isinf(a));
    if ~isempty(b);
        p01z(b,i)=min(p11z(c,i));
    end
end

[A_ga,B_ga,ADJga,ADJMga]=mos_cca([ga],mat,cs,cw,latx,laty,DIMS,PROP);
[A_gb,B_gb,ADJgb,ADJMgb]=mos_cca([gb],mat,cs,cw,latx,laty,DIMS,PROP);
if strcmp(cdib,'mixed_exp');
    [A_gc,B_gc,ADJgc,ADJMgc]=mos_cca([gc],mat,cs,cw,latx,laty,DIMS,PROP);
end
[A_01,B_01,ADJ01,ADJM01]=mos_cca([p01z],mat,cs,cw,latx,laty,DIMS,PROP);
[A_11,B_11,ADJ11,ADJM11]=mos_cca([p11z],mat,cs,cw,latx,laty,DIMS,PROP);

if nrun > 1
    [P01z,nb01]=search_mode(p01z,ADJM01);
    [P11z,nb11]=search_mode(p11z,ADJM11);
    [GA,nbga]=search_mode(ga,ADJMga);
    [GB,nbgb]=search_mode(gb,ADJMgb);
    if strcmp(cdib,'mixed_exp');
       [GC,nbgc]=search_mode(gc,ADJMgc);
       [GC]=scale_mean_var(GC,gc);
    end
    [GA]=scale_mean_var(GA,ga);
    [GB]=scale_mean_var(GB,gb);
    [P01z]=scale_mean_var(P01z,p01z);
    [P11z]=scale_mean_var(P11z,p11z);
else
    [P01z,nb01]=search_mode(p01z,ADJ01);
    [P11z,nb11]=search_mode(p11z,ADJ11);
    [GA,nbga]=search_mode(ga,ADJga);
    [GB,nbgb]=search_mode(gb,ADJgb);
    if strcmp(cdib,'mixed_exp');
       [GC,nbgc]=search_mode(gc,ADJgc);
       [GC]=scale_mean_var(GC,gc);
    end
    [GA]=scale_mean_var(GA,ga);
    [GB]=scale_mean_var(GB,gb);
    [P01z]=scale_mean_var(P01z,p01z);
    [P11z]=scale_mean_var(P11z,p11z);
end
    

for i=1:nc;
    a=find(GA(:,i) <= 0);
    b=find(GA(:,i) > 0);
    GA(a,i)=min(GA(b,i))*ones(size(a));
    a=find(GB(:,i) <= 0);
    b=find(GB(:,i) > 0);
    GB(a,i)=min(GB(b,i))*ones(size(a));
    if strcmp(cdib,'mixed_exp');
       a=find(GC(:,i) <= 0);
       b=find(GC(:,i) > 0);
       GC(a,i)=min(GC(b,i))*ones(size(a));
    end
end

P01=inverse_log_odd(P01z);
P11=inverse_log_odd(P11z);

for i=1:nyear;
    for j=1:nc;
        if strcmp(cdib,'gamma');
            [x_rand]=swg_gamma(-1,lseason,P01(i,j),P11(i,j),[GA(i,j) GB(i,j)],knn);
        elseif strcmp(cdib,'weibull');
            [x_rand]=swg_weibull(-1,lseason,P01(i,j),P11(i,j),[GA(i,j) GB(i,j)],knn);
        elseif strcmp(cdib,'mixed_exp');
            [x_rand]=swg_mixexp(-1,lseason,P01(i,j),P11(i,j),[GA(i,j) GB(i,j) GC(i,j)],knn);
        end
        X=x_rand';
        RAINMOD(((i-1)*lseason*knn)+1:i*lseason*knn,j)=X(:);
    end
end

        
