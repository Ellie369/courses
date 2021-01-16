function [NROW, NCOL, VALUE, NIR, NIC, FIR, FIC] = sparse_add_element(aij,i,j, NROW, NCOL, VALUE, NIR, NIC, FIR, FIC)
% %% for test this script
% clc; clear all;
% A=[1 3 4 8;2 1 2 3;4 3 5 8;9 2 7 4];
% % b=[1;1;1;1];
% [VALUE,FIC,FIR,NCOL,NROW,NIC,NIR]=Scheme0_storage(A)

%% add elements to VALUE, NROW, NCOL
n_VALUE=length(VALUE);
VALUE=[VALUE;aij];
NROW=[NROW;i];
NCOL=[NCOL;j];
num_aij=n_VALUE+1; % the number of aij in VALUE vector
index_R=FIR(i); 
index_C=FIC(j);

%% update FIR,NIR
if FIR(i)==0  %originally, all of the elements in ith row are zero, then aij would be the first element in ith row
    FIR(i)=num_aij; 
    NIR(num_aij)=0;
else %FIR(i)~=0
    if j<NCOL(index_R) % aij locates at the left side of FIR(i)
        FIR(i)=num_aij;
        NIR(num_aij)=index_R; 
    else %j>NCOL(index_R)
    while index_R~=0  % loop until 'index=the last element of ith row' 
        if NIR(index_R)==0 %index_R points to the last element of ith row
            NIR(index_R)=num_aij; 
            NIR(num_aij)=0;  
            break
        elseif (NCOL(index_R)<j)&&(j<NCOL(NIR(index_R))) %aij locates between two non-zero elements
            NIR(num_aij)=NIR(index_R);
            NIR(index_R)=num_aij;            
            break            
        end
            index_R=NIR(index_R);
    end                        
    end
end
%% update FIC NIC
if FIC(j)==0  %originally, all of the elements in jth row are zero, then aij would be the first element in jth column
    FIC(j)=num_aij; 
    NIC(num_aij)=0;
else %FIC(j)~=0
    if i<NROW(index_C) % aij locates above of FIC(i)
        FIC(j)=num_aij;
        NIC(num_aij)=index_C; 
    else %i>NROW(index_C)
    while index_C~=0  % loop until 'index=the last element of jth column' 
        if NIC(index_C)==0 %index_C points to the last element of jth column
            NIC(index_C)=num_aij; 
            NIC(num_aij)=0;  
            break
        elseif (NROW(index_C)<i)&&(i<NROW(NIC(index_C))) %aij locates between two non-zero elements
            NIC(num_aij)=NIC(index_C);
            NIC(index_C)=num_aij;            
            break            
        end
            index_C=NIC(index_C);
    end                        
    end
end    
end

            
            
            
            
            
      
    
        
 
    
    
    
    
    






