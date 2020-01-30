function [indice,l]=choose_date_daily(dep,year,month,day,hour);

% [indice]=choose_date_daily(depart,year,month,day,hour);
%
% This function gives indices of day in the daily series
% in any series beginning in date entered as 'dep.
%
% Input
% 'dep' : vector of four integers [year,month,day,hour] giving the first
% date of the file (hour should be 0 if the data are daily);
% 'year' : vector of two integers ([first year last year]) giving the years
% for which indices are needed
% 'month' : vector of two integers ([first month last month]) giving the
% months for which indices are needed
% 'day' : vector of two integers ([first day last day]) giving the
% day for which indices are needed
% 'hour' : one scalar (usually hour=0 if data are
% daily and could be 0,6,12 and 18 for other cases.
%
% ex : [indice]=choose_date([1957,9,1,0],[1961 1998],[6 10],[21 10],0); will gives
% the temporal indices of 21/6-10/10 of 1961-1998 in the ERA-40 dataset and
% [indice]=choose_date([1961,1,1,0],[1961 1998],[7 7],[1 31],12) will gives
% the temporal indices of 1-31/7 of 1961-1998 at 12 h).

N=datenum(dep(1),dep(2),dep(3),dep(4),0,0); % first date 
N=N-1;
indice=[];

for i=year(1):year(2);
    n_indice=length(indice);
    deb(i)=datenum(i,month(1),day(1),hour,0,0)-N;
    fin(i)=datenum(i,month(2),day(2),hour,0,0)-N;
    serie=[deb(i):1:fin(i)];
    l(i)=length(serie);
    indice(n_indice+1:n_indice+l(i))=serie(:);
end

