
 function [Scheme0_order]=Scheme0(A)
% clc; clear all;
% A=[1,1,0,1,0,0,0,1,0,0;1,1,0,0,0,0,1,0,0,1;0,0,1,1,1,0,1,1,0,1;1,0,1,1,1,0,0,0,0,0;0,0,1,1,1,1,1,0,0,1;0,0,0,0,1,1,1,0,0,0;0,1,1,0,1,1,1,1,1,0;1,0,1,0,0,0,1,1,0,0;0,0,0,0,0,0,1,0,1,0;0,1,1,0,1,0,0,0,0,1]
% spy(A)

%%
degrees=sum(A~=0)-1;
k=unique(degrees);
m=length(k);
index_re_degrees=[];
for i=1:m
    index_re_degrees=[index_re_degrees,find(degrees==k(i))]; 
end
Scheme0_order=index_re_degrees;
end
        
            
            
            
            
            
        
        
    
    

