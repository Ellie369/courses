
function [x]=LU_sparse(A_ordered,FIC,FIR,NCOL,NROW,NIC,NIR,b)
%% example for testing
% clc; clear all;
% A=[1 3 4 8;2 1 2 3;4 3 5 8;9 2 7 4];
% b=[1;1;1;1];
% [Scheme0_order]=Scheme0(A);
% [A_ordered,FIC,FIR,NCOL,NROW,NIC,NIR,b_ordered]=Sparse_storage(Scheme0_order, A,b)
%%
[QVALUE,QNROW, QNCOL, QNIR, QNIC, QFIR, QFIC] = Crout_Sparse(A_ordered,FIC,FIR,NCOL,NROW,NIC,NIR);
n=length(QFIR);
%% get vector y by forward substitution
y=zeros(n,1);
for k=1:n
    y_1=0;
    for j=1:k-1
        [Q_kj]=sparse_element(QVALUE,QFIC,QFIR,QNCOL,QNROW,QNIC,QNIR,k,j);        
        y_1= y_1+Q_kj*y(j);
    end
    [Q_kk]=sparse_element(QVALUE,QFIC,QFIR,QNCOL,QNROW,QNIC,QNIR,k,k);
    y(k)=(b(k)-y_1)/Q_kk;    
end
%% get the solution vector x by backward substitution
x=zeros(n,1);
for k=n:-1:1
    x_1=0;
    for j=k+1:n
         [Q_kj]=sparse_element(QVALUE,QFIC,QFIR,QNCOL,QNROW,QNIC,QNIR,k,j);
         x_1=x_1+Q_kj*x(j);
    end
    x(k)=y(k)-x_1;
end
end
