function[aij]=sparse_element(value,FIC,FIR,NCOL,NROW,NIC,NIR,i,j)
Index_value=FIR(i);
if Index_value==0  % FIR(i)=0,the i th row elements are all zeros and return aij=0; 
    aij = 0;
    return
else
    while Index_value~=0 %search all non-zero elements in ith row
        if NCOL(Index_value)==j
            aij=value(Index_value);
            return
        else
            Index_value=NIR(Index_value);
        end
        aij=0; % if not find a NCOL==j, then return 0
    end
end

    
            
        
    
    