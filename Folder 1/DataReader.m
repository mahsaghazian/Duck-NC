clear all
close all
%convert big data to profile number and depth and distance
for i=1:959
    filename=['File' num2str(i) '.csv'];
    Data = csvread(filename,0,1);
    p=Data(:,1);
    x=Data(:,7);
    d=Data(:,9);
    A=[p x d];
    filename1=['NewFile' num2str(i) '.xlsx' ];
    xlswrite(filename1,A); 
end
%seprate each profile data and intepolate to make each matrix the same
%length
for i=1:959
    filename=['NewFile' num2str(i) '.csv'];
    Data = csvread(filename,0,0);
    p=Data(:,1);
    if p==-91;
        Data(:,2)=x;
        Data(:,3)=d'
        
    end
    xq=(100:1500);
    dq=intrp1(x,d,xq);
     
end


