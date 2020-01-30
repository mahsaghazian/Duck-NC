function [dry_nb,dry_length,dry_sample,nb]=find_spell(X,thres,thres2,choice);

% [NB,MEAN,SAMPLE,ECH]=find_spell(X,thres,thres2,choice);
%
% Function to find spell in sequences of rainy days. 
%
% Input
% 'X' : vector of real number containing the rain for each day at one
% station. The column if any describes different years.
% 'thres' : any scalar to define the threshold from which a day is defined
% as wet (day > 'thres' is wet and a day receiving less or equal
% rainfall than 'thres' is dry). 
% 'thres2' : any scalar to efine the length of dry/wet spells for which the 
% user needs to know the umber 
% 'choice' : character string be set as 'dry' or 'wet'. If choice = 'wet',
% the results are about the wet sequences and if choice = 'dry', the
% results are about the dry sequences.
%
% Output
% 'NB' : matrix/vector of number of dry/wet spells
% 'MEAN' : matrix/vector of mean length of dry/wet spells
% 'SAMPLE' : matrix/vector of the whole distribution of dry/wet spells
% 'ECH' : matrix/vector of the number of wet/dry spells exceeding or equal 
% to 'thres2' consecutive days
%
% EXAMPLES
%
% A=[0 0 5.1 2 5 0 0 0 0 2 0 1 5 5 0 0 0 0 0 0 6 0.2 0.001 2 1 2 0 2 2 0]
%
% [a,b,c,d]=dryspell(A,0,6,'wet')
%
% a = 5
% b = 3
% c = 3 1 3 6 2
% d = 1
%
% [a,b,c,d]=dryspell(A,0,2,'dry')
% a = 6
% b = 2.5
% c = 2 4 1 6 1 1
% d = 3
%
% Vincent MORON
% september 2004
% modif : 4/5/2005 : Matrices are accepted as inputs (row describe day and
% column years)

[nday,nyear]=size(X);

if nargin==1; thres=0; end
Y=zeros(size(X));

if strcmp(choice,'wet');
   Y(find(X > thres))=ones(size(find(X > thres)));
elseif strcmp(choice,'dry');
   Y(find(X <= thres))=ones(size(find(X <= thres)));
end    

dry_sample=NaN*ones(size(Y));

for i=1:nyear;
    clear z zz zzz z2 dd
    dry_day=find(Y(:,i) == 1);
    nb_dry_day=length(dry_day);
    z=diff(dry_day);
    dry_nb(i)=nb_dry_day-length(find(z==1));
    if dry_nb(i)==1
        dry_sample(1,i)=(length(find(z==1))+1);
    elseif dry_nb(i) > 1
        zz=find(z > 1);
        zzz=diff(zz);
        z2=[zz(1);zzz(:)];
        dry_sample(1:dry_nb(i),i)=sort([z2;sum(nb_dry_day)-sum(z2)]);
    end
    dd=dry_sample(:,i);
    dd=dd(find(~isnan(dd)));
    dry_length(i)=sum(dd)/dry_nb(i);
    nb(i)=length(find(dd >= thres2));
end

dry_sample=dry_sample(1:max(dry_nb),:);
