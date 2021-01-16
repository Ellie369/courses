function [QVALUE,QNROW, QNCOL, QNIR, QNIC, QFIR, QFIC] = Crout_Sparse(A_ordered,FIC,FIR,NCOL,NROW,NIC,NIR)
%% only for testing this script
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
% [Scheme0_order]=Scheme0(A);
% [A_ordered,FIC,FIR,NCOL,NROW,NIC,NIR,b_ordered]=Sparse_storage(Scheme0_order, A,b);

% clc; clear all;
% A=[1 3 4 8;2 1 2 3;4 3 5 8;9 2 7 4];
% b=[1;1;1;1];
% [Scheme0_order]=Scheme0(A);
% [A_ordered,FIC,FIR,NCOL,NROW,NIC,NIR,b_ordered]=Sparse_storage(Scheme0_order, A,b)
%%
n=length(FIC); % get the dimention of A matrix  
QNROW=[];
QNCOL=[]; 
QVALUE=[];  
QNIR=[]; 
QNIC=[];  
QFIR=zeros(n,1); 
QFIC=zeros(n,1); 

% elements of jth column in Q matrix
for j=1:n  
    for k=j:n 
        [A_ordered_kj]=sparse_element(A_ordered,FIC,FIR,NCOL,NROW,NIC,NIR,k,j); % call an element in A_ordered 
        Qtemp_kj=0;        
        for i=1:j-1
            Q_ki=sparse_element(QVALUE,QFIC,QFIR,QNCOL,QNROW,QNIC,QNIR,k,i);
            Q_ij=sparse_element(QVALUE,QFIC,QFIR,QNCOL,QNROW,QNIC,QNIR,i,j);
            Qtemp_kj=Qtemp_kj+Q_ki*Q_ij;
        end        
        Q_kj=A_ordered_kj-Qtemp_kj;
        if Q_kj~=0
            [QNROW, QNCOL, QVALUE, QNIR, QNIC, QFIR, QFIC] = sparse_add_element(Q_kj, k, j, QNROW, QNCOL, QVALUE, QNIR, QNIC, QFIR, QFIC);
        end                      
    end
    % elements of jth row in Q matrix
    Q_jj=sparse_element(QVALUE,QFIC,QFIR,QNCOL,QNROW,QNIC,QNIR,j,j);
    if Q_jj~=0 
        for k=j+1:n
            Qtemp_jk=0;
            [A_ordered_jk]=sparse_element(A_ordered,FIC,FIR,NCOL,NROW,NIC,NIR,j,k);            
            for i=1:j-1
                Q_ji=sparse_element(QVALUE,QFIC,QFIR,QNCOL,QNROW,QNIC,QNIR,j,i);
                Q_ik=sparse_element(QVALUE,QFIC,QFIR,QNCOL,QNROW,QNIC,QNIR,i,k);
                Qtemp_jk=Qtemp_jk+Q_ji*Q_ik;
            end
            Q_jk=(A_ordered_jk-Qtemp_jk)/Q_jj;
            if Q_jk~=0
            [QNROW, QNCOL, QVALUE, QNIR, QNIC, QFIR, QFIC] = sparse_add_element(Q_jk, j,k, QNROW, QNCOL, QVALUE, QNIR, QNIC, QFIR, QFIC);
            end  
        end
    end
end

end
