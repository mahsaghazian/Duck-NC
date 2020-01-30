function [seas,occ,int,pw,pd,dry,wet]=subseas(obs,nseas,th);
        
% [seas,occ,int,pw,pd,dry,wet]=subseas(obs,nseas,th);
% 
% This function computes seeveral quantities characterizing the interannual
% variability of seasonal and subseasonal variations of daily rainfall.
%
% Input
% 'obs' : matrix of real number giving the daily rainfall at stations or
% grid-points. The rows describe time and the columns describe space.
% 'nseas' : integer scaler giving the length of each season (copied beneath
% ecah others in 'obs')
% 'th' : integer scalar to define wet day (if 'obs' > 'th')
%
% Output
% 'seas' : matrix of real number giving the seasonal amount at stations (as 
% the sum of daily rainfall)
% 'occ' : matrix of integer giving the number of rainy days (> 'th')
% 'int' : matrix of real number giving the daily mean intensity (for days >
% 'th')
% 'pw' : matrix of real number giving the persistence from wet to wet
% 'pd' : matrix of real number giving the persistence from dry to dry
% 'dry' : matrix of real number giving the mean length of dry spells (rain
% <= 'th')
% 'wet' : matrix of real number giving the mean length of wet spells (rain
% > 'th')
%
% Vincent MORON
% Dec 2004

[nt,nv]=size(obs);
[nyear]=nt/nseas;

for i=1:nyear
    for k=1:nv;
        O=obs(((i-1)*nseas)+1:i*nseas,k);
        seas(i,k)=sum(O);
        a=find(O > th);
        occ(i,k)=length(a);
        int(i,k)=sum(O(a))./occ(i,k);
        pw(i,k)=length(find(O(1:nseas-1) > th & O(2:nseas) > th))./length(find(O(1:nseas-1) > th));
        pd(i,k)=length(find(O(1:nseas-1) <= th & O(2:nseas) <= th))./length(find(O(1:nseas-1) <= th));
        [dd,dry(i,k),dry2,dd3]=find_spell(O,th,1,'dry');
        [dd,wet(i,k),dry2,dd3]=find_spell(O,th,1,'wet');
    end
end
