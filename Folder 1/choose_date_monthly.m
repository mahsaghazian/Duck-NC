function [indice]=choose_date_monthly(year1,year2,nmon,choice);

% [indice]=choose_date_monthly(year1,year2,nmon);
%
% this function creates temporal indices from monthly time scales data (for
% daily time scale, use choose_data_daily.m).
%
% Input
% 'year1' : a vector of two integers giving the first and last year of the
% total period
% 'year2' : a vector of two integers giving the first and last year of the
% period to be selected
% 'nmon' : a vector of integer giving the month to be selected (if you want
% select a season across a year like December-February, nmon = [12 1 2] and
% NOT [1 2 12].
% 'choice' : string character to define the fact that the output are
% complete seasons or not in case of a season across a year. If choice =
% 'no' and nmon=[12 1 2], the January-February of the firt year and the 
% December of the last year are removed and kept if choice='yes'.
% 
% Output
% 'indice' : a column vector of integer giving the indices choosen. If the
% months selected are across the year like DJF and choice='no', the first indice is
% december of the first year, the second indice is january of the second
% year, so that no indice are given for January-February of the first year
% and December of the last year
%
% Vincent Moron
% Nov 2005

nyear=(year1(2)-year1(1))+1;
d1year=year2(1)-year1(1)+1;
d2year=year2(2)-year1(1)+1;
len=nyear*12;
len2=((year2(2)-year2(1))+1)*length(nmon);
a=reshape([1:nyear*12],12,nyear);
nmon2=sort(nmon);
if (~isempty(find(nmon==12)) & ~isempty(find(nmon==1)));
    c=find(diff(nmon2) > 1)+1;
    d=len2-find(diff(nmon) < 0);
else
    c=1;
    d=len2;
end
indice=a(nmon2,d1year:d2year);
indice=indice(:);
if strcmp(choice,'no');
    indice=indice(c:d);
end
    