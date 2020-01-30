clear all
close all
%Date
for i=1:912
    filename=['File' num2str(i) '.csv'];
    t=csvread(filename,0,11);
    t1=t(1,1);
    TT(i,1)=t1;
end

figure
n=0;
m=1;
for i=1:912
    j=1;
    filename=['Profile' num2str(i) ',' num2str(j) '.csv'];
    if exist(filename, 'file')
     n=n+1;
    x=csvread(filename,0,0);
    d=csvread(filename,0,1);
    dx=d(:,1);
    
    Tm(m)=TT(i);
    Tq(m)=i;
    Cq=[Tm ;Tq];

   x1=x(:,1);
   m1=700;
    xq=(100:m1);
    [x1, index] = unique(x1); 
dq = interp1(x1, dx(index), xq);
    if isnan(dq)==0
     m=m+1;
    end
    h(n,:)=dq;
    hq=denan(h);
    
    end
end

   
save('Time.mat', 'Cq')
example = matfile('Profile.no.mat');
q = example.U;
q=q(1:408);
    
    for i=1:408
    j=1;
   
    k=q(i);
    filename2=['Profile' num2str(k) ',' num2str(j) '.csv'];
    if exist(filename2, 'file')
    
    x=csvread(filename2,0,0);
    d=csvread(filename2,0,1);
    dx=d(:,1);
   x1=x(:,1);
   
    xq=(100:m1);
    [x1, index] = unique(x1); 
dq = interp1(x1, dx(index), xq); 
    H(i,:)=dq;
    
    hold on
    plot(xq,dq);
    xlabel('Distance(m)','FontWeight','bold','fontsize',16)
    ylabel('Elevation(m)','FontWeight','bold','fontsize',16)
    end
end
  
   
l1=size(H);



figure;
surface(xq,q,H,'EdgeColor','none','LineStyle','none');

  