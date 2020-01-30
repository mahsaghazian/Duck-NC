function [lepsp,leps,lepspotts,potts,SK]=lepscont2(yobs,yestim);

% [lepsp,leps,lepspotts,potts,,SK]=lepscont2(yobs,yestim);
% 
% This function computes linear error in painting space (LEPS) for
% continuous forecast.
%
% Input
% 'yobs' : a vector of real number giving the observations
% 'yestim' : a vector of real number giving the forecast
% 
% Output
% 'lepsp' : integer giving the mean LEPS in %
% 'leps' : vector of real number giving the LEPS for each time unit.
% 'lepspotts : integer giving the mean revised LEPS 
% 'potts' : vector of real number giving the revised LEPS for each time
% unit
% 'SK' : integer giving the skill (between -100 and 100 from 'potts' ; a
% value of 100% means a perfect forecast, that is the rank of 'yobs' is
% perfectly match by the one of 'yestim'.
%
% Nathalie Philippon & Pascal Oettli
% Fev 2003

yobs=yobs(:);
yestim=yestim(:);
MU = mean(yobs); SIGMA = std(yobs);
Pv = normcdf(yobs,MU,SIGMA);
Pf = normcdf(yestim,MU,SIGMA);

a = abs(Pf-Pv);
S = 1-a;
Cs = 1.5 .* (Pf-(Pf.^2)+0.5) .* (Pv-(Pv.^2)+0.5);
leps = S-Cs;
Sp = 1-Cs;
lepsp = (sum(leps)./sum(Sp)).*100;

potts=3.*(1-abs(Pf-Pv)+(Pf).^2-Pf+(Pv).^2-Pv)-1;
lepspotts=mean(potts);

%%% si potts >0
potts2=3.*(1-abs(Pv-Pv)+(Pv).^2-Pv+(Pv).^2-Pv)-1;

%%% si potts <0
potts3=3.*(1-abs(1-Pv)+1.^2-1+(Pv).^2-Pv)-1;
potts4=3.*(1-abs(0-Pv)+0.^2-0+(Pv).^2-Pv)-1;
potts34=min([potts3 potts4],[],2);

num=sum(100.*potts);

if num>0,
SK=(num)/(sum(potts2));
else,
SK=(num)/(sum(abs(potts34)));
end


