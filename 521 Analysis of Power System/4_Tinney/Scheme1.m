
 function [Scheme1_order]=Scheme1(A)
% clc; clear all;
% 10 nodes example
% A=[1,1,0,1,0,0,0,1,0,0;1,1,0,0,0,0,1,0,0,1;0,0,1,1,1,0,1,1,0,1;1,0,1,1,1,0,0,0,0,0;0,0,1,1,1,1,1,0,0,1;0,0,0,0,1,1,1,0,0,0;0,1,1,0,1,1,1,1,1,0;1,0,1,0,0,0,1,1,0,0;0,0,0,0,0,0,1,0,1,0;0,1,1,0,1,0,0,0,0,1]
degrees=sum(A~=0)-1;   %set the start value for degrees
Scheme1_order=[];
% A_ordered=A;
while max(degrees)>0     
    index_re_degrees=[];     
    k=unique(degrees);
    k_index=find(k>0);
    m=length(k_index);
    for i=1:m
        index_re_degrees=[index_re_degrees,find(degrees==k(k_index(i)))];
    end
    ii=index_re_degrees(1);
    Scheme1_order=[Scheme1_order,ii];    
% find non-zero elements in ii row of A matrix which has the least degree
    A1_nnz=find(A(ii,:)~=0);
    if length(A1_nnz)>2 
        for i=1:length(A1_nnz)-1    
            for j=i+1:length(A1_nnz)
                if (A1_nnz(i)~=ii)&&(A1_nnz(j)~=ii)&&A(A1_nnz(i), A1_nnz(j))==0  %build connection for other nodes
                    A(A1_nnz(i), A1_nnz(j))=1;
                    A(A1_nnz(j), A1_nnz(i))=1;
                end
            end
        end
    end
    A(ii,:)=0; %delete the node with the least degree
    A(:,ii)=0;
    degrees=sum(A~=0)-1;
end 
    Scheme1_order=[Scheme1_order,index_re_degrees(end)];  
end


        
            
            
            
            
            
        
        
    
    

