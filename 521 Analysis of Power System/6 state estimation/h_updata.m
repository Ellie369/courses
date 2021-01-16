function [h]=h_updata(V_mag,V_Delta,Line_data,No_of_Buses,No_of_Lines,Theta,Y_mag)
                   
% P_Calculation and Q_Calculation update 
[P_cal,Q_cal]=cal_PQ(V_mag,Y_mag,Theta,V_Delta,No_of_Buses); 
 
% Pij;Qij;Pji;Qji update 
Pij = zeros(No_of_Lines,1);  Qij = zeros(No_of_Lines,1);  Pji = zeros(No_of_Lines,1);  Qji = zeros(No_of_Lines,1);
for k = 1:No_of_Lines
    i = Line_data(k,2); j = Line_data(k,3);
    Pij(k) = V_mag(i)*V_mag(j)*Y_mag(i,j)*cos(V_Delta(i)-V_Delta(j)-Theta(i,j)) - V_mag(i)^2*Y_mag(i,j)*cos(Theta(i,j));
    Qij(k) = V_mag(i)*V_mag(j)*Y_mag(i,j)*sin(V_Delta(i)-V_Delta(j)-Theta(i,j)) + V_mag(i)^2*Y_mag(i,j)*sin(Theta(i,j));
    Pji(k) = V_mag(j)*V_mag(i)*Y_mag(j,i)*cos(V_Delta(j)-V_Delta(i)-Theta(j,i)) - V_mag(j)^2*Y_mag(j,i)*cos(Theta(j,i));
    Qji(k) = V_mag(j)*V_mag(i)*Y_mag(j,i)*sin(V_Delta(j)-V_Delta(i)-Theta(j,i)) + V_mag(j)^2*Y_mag(j,i)*sin(Theta(j,i));
end
h = [V_mag;P_cal;Q_cal;Pij;Pji;Qij;Qji];