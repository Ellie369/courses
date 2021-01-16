
 function [A_ordered,FIC,FIR,NCOL,NROW,NIC,NIR,b_ordered]=Sparse_storage(index_re_degrees, A,b)
% clc; clear all;
% A=[11 12 13 14 15;21 22 0 0 0;31 0 33 0 0;41 0 0 44 0;51 0 0 0 55];
% b=[1;2;0;4;5];
% A=[1 3 4 8;2 1 2 3;4 3 5 8;9 2 7 4];
% b=[1;1;1;1];
% [index_re_degrees]=Scheme0(A);

%% storage the matrix using sparse technique
A_ordered=[]; % the non zero elements in A matrix
% b_ordered=[];
n=length(A);
NROW=[];
NCOL=[];
for i=1:n 
    ii=index_re_degrees(i);
    for j=1:n
        jj=index_re_degrees(j);
        if A(ii,jj)~=0
            A_ordered=[A_ordered; A(ii,jj)];
            NROW=[NROW;i];
            NCOL=[NCOL;j];                     
        end
    end
end
nnz=length(A_ordered);
NIR=zeros(nnz,1);
NIC=zeros(nnz,1);
for i=1:nnz-1
    if NROW(i+1)==NROW(i)
       NIR(i)=i+1;         
    end   
end

for i=1:nnz-1
    for j=i+1:nnz
        if NCOL(j)==NCOL(i)
            NIC(i)=j;
            break
        end        
    end
end

FIR=zeros(n,1);
FIC=zeros(n,1);
for i=1:n
    for j=1:nnz   
        if NROW(j)==i
        FIR(i)=j;
        break
        end        
    end
end
for i=1:n
    for j=1:nnz   
        if NCOL(j)==i
        FIC(i)=j;
        break
        end        
    end
end
%% order b vector
b_ordered=[];
for i=1:n 
    ii=index_re_degrees(i);
    b_ordered=[b_ordered; b(ii)]; 
end
end
        
            
            
            
            
            
        
        
    
    

