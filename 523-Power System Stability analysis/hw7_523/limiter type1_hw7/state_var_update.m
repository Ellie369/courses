% x(1)...x(11) --------- voltage magnitudes
% x(12)...x(22) -------- voltage angles

function f= state_var_update(x,No_of_Buses,Y_ang,Y_mag,PL0, M0,G0,QL0,H0,B0,P_gen,Q_gen)
%% make the slack bus voltage and angle fixed

x(1) = 1.03;  
x(1+No_of_Buses) = 0;

%% Calculate Pi and Qi
P_net=zeros(No_of_Buses,1);
Q_net=zeros(No_of_Buses,1);
for i=1:No_of_Buses
    for j=1:No_of_Buses
    P_net(i)=P_net(i)+x(i)*x(j)*Y_mag(i,j)*cos(x(No_of_Buses+i)-x(No_of_Buses+j)-Y_ang(i,j));
    Q_net(i)=Q_net(i)+x(i)*x(j)*Y_mag(i,j)*sin(x(No_of_Buses+i)-x(No_of_Buses+j)-Y_ang(i,j));
    end 
end

%% load
P_load =PL0+(M0).*x(1:No_of_Buses) + (G0).* (x(1:No_of_Buses).^2);
Q_load =QL0+(H0).*x(1:No_of_Buses) + (B0).* (x(1:No_of_Buses).^2);

%% power flow equations
% f1-active power equation; f2-reactive power equation
f1 = (P_gen - P_load - P_net);
f2 = (Q_gen - Q_load- Q_net);
f=[f1;f2];
end