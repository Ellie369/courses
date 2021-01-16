function H=H_Matrix(V_mag,V_Delta,Line_data,No_of_Buses,No_of_Lines,Theta,Y_mag,B,G)
% H1 and H2
H1 = zeros(No_of_Buses, No_of_Buses-1);  H2 = eye(No_of_Buses);

% H3 to H6

H3 = zeros(No_of_Buses, No_of_Buses-1);  H4 = zeros(No_of_Buses);
H5 = zeros(No_of_Buses, No_of_Buses-1);  H6 = zeros(No_of_Buses);

for i = 1:No_of_Buses
            for j = 2:No_of_Buses
                if i~=j
                    H3(i,j-1) = -abs(V_mag(i)*V_mag(j)*Y_mag(i,j))*sin(Theta(i,j)+V_Delta(j)-V_Delta(i));
                    H5(i,j-1) = -abs(V_mag(i)*V_mag(j)*Y_mag(i,j))*cos(Theta(i,j)+V_Delta(j)-V_Delta(i));
                else
                    for n = [1:i-1, i+1:No_of_Buses]
                        H3(i,j-1) = H3(i,j-1) + abs(V_mag(i)*V_mag(n)*Y_mag(i,n))*sin(Theta(i,n)+V_Delta(n)-V_Delta(i));
                        H5(i,j-1) = H5(i,j-1) + abs(V_mag(i)*V_mag(n)*Y_mag(i,n))*cos(Theta(i,n)+V_Delta(n)-V_Delta(i));
                    end
                end
            end
        end
        %  H4 H6
        for i = 1:No_of_Buses
            for j = 1:No_of_Buses
                if i~=j
                    H4(i,j) = abs(V_mag(i)*Y_mag(i,j))*cos(Theta(i,j)+V_Delta(j)-V_Delta(i));
                    H6(i,j) = -abs(V_mag(i)*Y_mag(i,j))*sin(Theta(i,j)+V_Delta(j)-V_Delta(i));
                else
                    for n = [1:i-1, i+1:No_of_Buses]
                        H4(i,j) = H4(i,j) + abs(V_mag(n)*Y_mag(i,n))*cos(Theta(i,n)+V_Delta(n)-V_Delta(i));
                        H6(i,j) = H6(i,j) + abs(V_mag(n)*Y_mag(i,n))*sin(Theta(i,n)+V_Delta(n)-V_Delta(i));
                    end
                    H4(i,j) = H4(i,j) + 2*V_mag(i)*G(i,j);
                    H6(i,j) = -H6(i,j) - 2*V_mag(i)*B(i,j);
                end
            end
        end
       

% H7 to H14

H7 = zeros(No_of_Lines, No_of_Buses);  H8 = zeros(No_of_Lines, No_of_Buses);
H9 = zeros(No_of_Lines, No_of_Buses);  H10 = zeros(No_of_Lines, No_of_Buses);
H11 = zeros(No_of_Lines, No_of_Buses);  H12 = zeros(No_of_Lines, No_of_Buses);
H13 = zeros(No_of_Lines, No_of_Buses);  H14 = zeros(No_of_Lines, No_of_Buses);

for k = 1:No_of_Lines
            i = Line_data(k,2); j = Line_data(k,3);
            H7(k,j) = -abs(V_mag(i)*V_mag(j)*Y_mag(i,j))*sin(Theta(i,j)+V_Delta(j)-V_Delta(i));
            H7(k,i) = abs(V_mag(i)*V_mag(j)*Y_mag(i,j))*sin(Theta(i,j)+V_Delta(j)-V_Delta(i));
            H8(k,j) = abs(V_mag(i)*Y_mag(i,j))*cos(Theta(i,j)+V_Delta(j)-V_Delta(i));
            H8(k,i) = -2*abs(V_mag(i))*G(i,j) + abs(V_mag(j)*Y_mag(i,j))*cos(Theta(i,j)+V_Delta(j)-V_Delta(i));
            H9(k,j) = abs(V_mag(j)*V_mag(i)*Y_mag(j,i))*sin(Theta(j,i)+V_Delta(i)-V_Delta(j));
            H9(k,i) = -abs(V_mag(j)*V_mag(i)*Y_mag(j,i))*sin(Theta(j,i)+V_Delta(i)-V_Delta(j));
            H10(k,j) = -2*abs(V_mag(j))*G(j,i) + abs(V_mag(i)*Y_mag(j,i))*cos(Theta(j,i)+V_Delta(i)-V_Delta(j));
            H10(k,i) = abs(V_mag(j)*Y_mag(j,i))*cos(Theta(j,i)+V_Delta(i)-V_Delta(j));
            H11(k,j) = -abs(V_mag(i)*V_mag(j)*Y_mag(i,j))*cos(Theta(i,j)+V_Delta(j)-V_Delta(i));
            H11(k,i) = abs(V_mag(i)*V_mag(j)*Y_mag(i,j))*cos(Theta(i,j)+V_Delta(j)-V_Delta(i));
            H12(k,j) = -abs(V_mag(i)*Y_mag(i,j))*sin(Theta(i,j)+V_Delta(j)-V_Delta(i));
            H12(k,i) = -2*abs(V_mag(i))*(-B(i,j)) - abs(V_mag(j)*Y_mag(i,j))*sin(Theta(i,j)+V_Delta(j)-V_Delta(i));
            H13(k,j) = abs(V_mag(j)*V_mag(i)*Y_mag(j,i))*cos(Theta(j,i)+V_Delta(i)-V_Delta(j));
            H13(k,i) = -abs(V_mag(j)*V_mag(i)*Y_mag(j,i))*cos(Theta(j,i)+V_Delta(i)-V_Delta(j));
            H14(k,j) = -2*abs(V_mag(j))*(-B(j,i)) - abs(V_mag(i)*Y_mag(j,i))*sin(Theta(j,i)+V_Delta(i)-V_Delta(j));
            H14(k,i) = -abs(V_mag(j)*Y_mag(j,i))*sin(Theta(j,i)+V_Delta(i)-V_Delta(j));
        end
        H7(:,1) = []; H9(:,1) = [];  H11(:,1) = [];  H13(:,1) = [];
 %  H matrix
 H = [H1 H2; H3 H4; H5 H6; H7 H8; H9 H10; H11 H12; H13 H14];
end