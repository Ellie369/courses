function [Q, alpha,beta]=Qmatrix_alpha_beta(A)
% % only for testing this script
% clc; clear all;
% % 10 nodes example
% A=[11 12 0 14 0 0 0 18 0  0; %1
%    21 22 0 0 0 0 27 0 0 210; %2
%    0 0 33 34 35 0 37 38 0 310; %3 
%    41 0 43 44 45 0 0 0 0  0; %4
%    0 0 53 54 55 56 57 0 0 510; %5
%    0 0 0 0 65 66 67 0 0  0; %6
%    0 1 2 0 75 76 77 78 79  0; %7
%    1 0 83 0 0 0 87 88 0  0; %8
%    0 0 0 0 0 0 97 0 99  0; %9
%    0 102 103 0 105 0 0 0 0 1010; %10
%    ];
% b=[1 2 3 4 5 6 7 8 9 10]';
%%
n=length(A); % get the dimension of A matrix    
Q=zeros(n);
Q_1=zeros(n); 
alpha_col=zeros(n,1);
alpha_row=zeros(n,1);
alpha=0;
beta=0;
for j=1:n
% elements of jth column in Q matrix
    for k=j:n 
        for i=1:j-1
            Q_1(k,j)=Q_1(k,j)+Q(k,i)*Q(i,j);
        end
        Q(k,j)=A(k,j)-Q_1(k,j);
        if Q(k,j)~=0
            alpha_col(j)=alpha_col(j)+1;
            beta=beta+1;
        end
     end
% elements of jth row in Q matrix
if Q(j,j)~=0     
for k=j+1:n    
    for i=1:j-1
        Q_1(j,k)=Q_1(j,k)+Q(j,i)*Q(i,k);
    end
    Q(j,k)=(A(j,k)-Q_1(j,k))/Q(j,j);
    if Q(j,k)~=0
        alpha_row(j)=alpha_row(j)+1;    
        beta=beta+1;
    end
end
end
alpha=alpha+alpha_col(j)*alpha_row(j);
end
end





