%% read data
clc;
close all
clear all;
%%%%%%%%%%%%%%
%% Part1 Economic Dispatch without power flow
%%%%%%%%%%%%%%
Ptotal=259;
A=[0.008,0,-1;0 0.0096 -1;1 1 0]; 
b = [-8;-6.4;Ptotal];
[x]=LU_factor_PQ(A,b);
P1 = x(1);  
P2 = x(2); 
f1=0.004 * P1^2 + 8 * P1;
f2=0.0048 * P2^2 + 6.4 * P2;
f=f1+f2;
Economic_dispatch_1=[P1;P2;f];
fprintf('Ecnomic dispatch for %4.3f MW load:\n',Ptotal);
fprintf('The output of P1 and P2 are %4.3f MW, %4.3f MW and the total cost is %4.3f\n',P1,P2,f);

%%%%%%%%%%%%%%%
%%Part 2 OPF with u=P2, and P12<=5MW
%%%%%%%%%%%%%%%

 %% read grid data from file
file_name='ieee14cdf.txt';
[S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data]=read_data(file_name);
[Y_mat,Theta,Y_mag,B,G]=y_bus(Bus_data,Line_data,No_of_Buses,No_of_Lines); %form Y_matrix
P2=P2/100;
P12_max=5; % P12<=5MW;
% P12_max=1e5;
gradient_L_u =1; % initia value for gradient_L_u
stepsize = 1;  %step size
iter=0;

while abs(gradient_L_u) >1e-3
    [PQ,nPQ,P_gen_cal,J,V_mag,V_Delta]= NRpowerflow(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data,P2);
    [gradient_L_u,P12]=gradient_method(P1,P2,P12_max,J,V_mag,V_Delta,Theta,Y_mag,PQ,nPQ);
    P2 = P2 - stepsize * gradient_L_u/100;
    P1 = P_gen_cal(1)*100;
    iter=iter+1;
end

P2=P2*100;
f1=0.004 * P1^2 + 8 * P1;
f2=0.0048 * P2^2 + 6.4 * P2;
f=f1+f2;
fprintf('Result for optimal power flow:\n');
fprintf('The output of P1 is %4.3f MW, P2 is %4.3f MW,P12 flow equals to %4.3f MW and the total cost is %4.3f\n',P1,P2,P12,f);
%%%%%%%%%%%%%%%

