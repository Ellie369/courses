function[B1,B2]=B_mat(B,nPQ,PQ,No_of_Buses)
B1=-B(2:No_of_Buses,2:No_of_Buses);

for i=1:nPQ
    for j=1:nPQ
        B2(i,j)=-B(PQ(i),PQ(j));
    end
end
end



