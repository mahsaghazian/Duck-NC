function [adjusted,nomb]=search_mode(obs,sim);

% [adjusted,nb]=search_mode(obs,sim);
% 
% THis function find the "best" number of CCA modes to be included in the 
% reconstruction of a matrix 'obs' by a set of modes in 'sim'. It is used 
% with mos_cca.m to select the number of CCA modes to be retained in the
% reconstruction. The "best" number 'nb' is the one for which the median of
% the correlations between 'obs' and 'sim' is highest.
%
% Input
% 'obs' : a matrix of real number with rows describing time and columns
% describing space.
% 'sim' : a matrix of real number with the reconstruction of 'obs' with the 
% CCA modes. the first columns are the reconstruction of 'obs' with the 
% leading CCA mode, then the next ones are the reconstruction with the 
% second one and so on
%
% Output
% 'adjusted' : a matrix of real number of the sum of the reconstructions
% (same size as 'obs');
% 'nb' : scalr integer giving the number of modes used in the
% reconstruction
%
% Vincent Moron
% Nov 2005

[nr,nc]=size(obs);
[nr,nc2]=size(sim);
nb=nc2/nc;
R(:,1)=diag((1./(nr-1))*stan(obs,'s')'*stan(sim(:,1:nc),'s'));

if nb > 1
    for i=2:nb;
        S=reshape(sim(:,1:i*nc),nc*nr,i);
        S=sum(S');
        S=reshape(S,nr,nc);
        R(:,i)=diag((1./(nr-1))*stan(obs,'s')'*stan(S,'s'));
    end
end

nomb=find(median(R)==max(median(R)));
adjusted=reshape(sim(:,1:nomb*nc),nr*nc,nomb);
if nomb > 1;
    adjusted=sum(adjusted');
end
adjusted=reshape(adjusted,nr,nc);
