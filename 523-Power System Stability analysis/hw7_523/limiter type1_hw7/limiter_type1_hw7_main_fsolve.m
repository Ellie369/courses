%% homework 7
%generator parameters, changing rating from 900MVA to 100MVA 
%p---prime; pp---double prime
clc
clear all
 format shortEng
 format compact
 %% Type1 model fault
fault_time =2;    %%t=t+1   %% This is the time when the fault occurs
clearing_time = 3/60;        %% This is the time of clearance after the fault has occurred
endtime = 21;    %% t=t+1
h = 1/60;   %% Step size for simulation

X_d(2:4,1)=1.8*100/900;
X_q(2:4,1)=1.7*100/900;
X_l(2:4,1)=0.2*100/900;
X_dp(2:4,1)=0.3*100/900;
X_qp(2:4,1)=0.55*100/900;
X_dpp(2:4,1)=0.25*100/900;
X_qpp(2:4,1)=0.25*100/900;
R_a(2:4,1)=0; 
T_d0p(2:4,1)=8.0;
T_q0p(2:4,1)=0.4;
T_d0pp(2:4,1)=0.03;
T_q0pp(2:4,1)=0.05;
H(1)=6.5*900/100;
H(2)=6.5*900/100;
H(3)=6.175*900/100;
H(4)=6.175*900/100;
K_D(2:4,1)=2*900/100;

X=(X_d+X_q)/2;
X_p=(X_dp+X_qp)/2;
ws=2*pi*60;
no_of_states = 6; 
%% Exciter parameters
KA=50;
TA=0.01;
VRmin=-4;
VRmax=4;
Efdmin=0;
Efdmax=2.0;

%% Governer parameters
Tsg=100;
Ksg=1;
Psgmin=0;
Psgmax=1;
R=0.0056; %convert to Base =100MVA
%% read grid data from file
file_name='b_kundur_system.txt';
[S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data]=read_data(file_name);
 [PQ,nPQ,PV,nPV,Y_ang,Y_mag,B,G,V_mag,V_Delta,P_gen,Q_gen, V_result,P_load, Q_load]= NR_power_flow_type1(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data);

%% solve equilibrium points for dynamic equations
V_bus0=V_mag;
Delta0=V_Delta;
for genbus=1:4 %generator buses, including slack bus
    x0=[0;1;1;0.5;1;1;1.7;7;1;7;7]; %
    options=optimoptions('fsolve','algorithm','levenberg-marquardt','display','off');
    [solution,~,exitflag]  = fsolve(@(x) type1_equilibrium_points(x,genbus,ws,K_D,X_d,X_dp,X_q,X_qp,KA,P_gen, Q_gen, V_bus0,Delta0,R),x0,options);
    theta_0(genbus,1)=solution(1);
    omega_0(genbus,1)=solution(2);
    E_qp0(genbus,1)=solution(3);
    E_dp0(genbus,1)=solution(4);
    E_fd0(genbus,1)=solution(7);
    Pm_0(genbus,1)=solution(8);
    V_ref(genbus,1)=solution(9);
    Pc_0(genbus,1)=solution(10);
    Pe_0(genbus,1)=solution(11);
end
%% load parameters
%% Calculate Pload
%%consider load as ZIP model
% initialize ZIP lodad model
P_p=0.5; P_q=0.5; %constant power percent
I_p=0; I_q=0;   %constant current percent
Z_p=0.5; Z_q=0.5; % constant impedance percent

PL0=P_p* P_load;
M0=(I_p* P_load)./V_mag;
G0=(Z_p*P_load)./(V_mag.^2);
QL0=P_q* Q_load;
H0=(I_q* Q_load)./V_mag; 
B0=(Z_q*Q_load)./(V_mag.^2);


%% pre fault, from t=0 to t=fault-
iter=1;
for iteration=1:h:fault_time
    % update state variables
    x0= [V_mag; V_Delta]; %% Initialize the starting point
     options=optimoptions('fsolve','display','off');
    [solution,~,exitflag] = fsolve(@(x) state_var_update(x,No_of_Buses,Y_ang,Y_mag, PL0, M0,G0,QL0,H0,B0,P_gen,Q_gen),x0,options);
    
    V_mag = solution(1:No_of_Buses); % obtain voltages from fsolve solution
    V_Delta = solution(No_of_Buses+1:end); % obtain angles from fsolve solution
% update dynamic state variables
 [theta_0,omega_0,E_dp0, E_qp0,E_fd0, Pm_0,P_gen,Q_gen,V_mag]=type1_euler(Pm_0, P_gen,Q_gen, Pc_0,H,V_mag,V_Delta, theta_0,omega_0,ws,K_D,h,...
        E_qp0, E_dp0, E_fd0, V_ref, X_d, X_q, X_dp, X_qp, T_d0p,T_q0p,KA,TA,Tsg,Ksg,R,VRmax,VRmin,Efdmax,Efdmin,Psgmax,Psgmin);
    time(iter,1)=iteration-1;
    theta(:,iter)=theta_0;
    w(:,iter)=omega_0;
    Edp(:,iter)= E_dp0;
    Eqp(:,iter)= E_qp0;
    Efd(:,iter)= E_fd0;
    Pm(:,iter)= Pm_0;
    iter=iter+1;
end
iter=iter-1; %% To delete the extra time step added at the end
%     plot(time,rad2deg(theta(2,:)),'linewidth',1.5)
 plot(time,Pm)

%% on fault
% revise y matrix, and form y_gen
file_name_2='onfault_b_kundur_system.txt';
[S_Base_2,No_of_Buses_2,No_of_Lines_2,Bus_data_2,Line_data_2]=read_data(file_name_2);
[Y_mat_f,Y_ang_f,Y_mag_f,B_f,G_f]=y_bus(Bus_data_2,Line_data_2,No_of_Buses_2,No_of_Lines_2);
% V_mag=[V_mag;1]; 
% V_Delta=[V_Delta;0];

for iteration=fault_time:h:(fault_time+clearing_time)
    % update state variables
    x0= [V_mag; V_Delta]; %% Initialize the starting point
    options=optimoptions('fsolve','display','off');
    [solution,~,exitflag] = fsolve(@(x) state_var_update(x,No_of_Buses_2,Y_ang_f,Y_mag_f,PL0, M0,G0,QL0,H0,B0,P_gen,Q_gen),x0,options);
    
    V_mag = solution(1:No_of_Buses_2); % obtain voltages from fsolve solution
    V_Delta = solution(No_of_Buses_2+1:end); % obtain angles from fsolve solution
% update dynamic state variables
 [theta_0,omega_0,E_dp0, E_qp0,E_fd0, Pm_0,P_gen,Q_gen,V_mag]=type1_euler(Pm_0, P_gen,Q_gen, Pc_0,H,V_mag,V_Delta, theta_0,omega_0,ws,K_D,h,...
        E_qp0, E_dp0, E_fd0, V_ref, X_d, X_q, X_dp, X_qp, T_d0p,T_q0p,KA,TA,Tsg,Ksg,R,VRmax,VRmin,Efdmax,Efdmin,Psgmax,Psgmin);
   
    time(iter,1)=iteration-1;
    theta(:,iter)=theta_0;
    w(:,iter)=omega_0;
    Edp(:,iter)= E_dp0;
    Eqp(:,iter)= E_qp0;
    Efd(:,iter)= E_fd0;
    Pm(:,iter)= Pm_0;
    iter=iter+1;
end
iter=iter-1; %% To delete the extra time step added at the end
%     plot(time,rad2deg(theta(2,:)),'linewidth',1.5)
 plot(time,Pm)
%% post fault analysis
% revise y matrix
file_name_3='postfault_b_kundur_system.txt';
[S_Base_3,No_of_Buses_3,No_of_Lines_3,Bus_data_3,Line_data_3]=read_data(file_name_3);
[Y_mat_pf,Y_ang_pf,Y_mag_pf,B_pf,G_pf]=y_bus(Bus_data_3,Line_data_3,No_of_Buses_3,No_of_Lines_3);

for iteration= (fault_time+clearing_time):h:endtime
    % update state variables
    x0= [V_mag; V_Delta]; %% Initialize the starting point
    options=optimoptions('fsolve','display','off');
    [solution,~,exitflag] = fsolve(@(x) state_var_update(x,No_of_Buses_3,Y_ang_pf,Y_mag_pf,PL0, M0,G0,QL0,H0,B0,P_gen,Q_gen),x0,options);
    
    V_mag = solution(1:No_of_Buses); % obtain voltages from fsolve solution
    V_Delta = solution(No_of_Buses+1:end); % obtain angles from fsolve solution
% update dynamic state variables
  [theta_0,omega_0,E_dp0, E_qp0,E_fd0, Pm_0,P_gen,Q_gen,V_mag]=type1_euler(Pm_0, P_gen,Q_gen, Pc_0,H,V_mag,V_Delta, theta_0,omega_0,ws,K_D,h,...
        E_qp0, E_dp0, E_fd0, V_ref, X_d, X_q, X_dp, X_qp, T_d0p,T_q0p,KA,TA,Tsg,Ksg,R,VRmax,VRmin,Efdmax,Efdmin,Psgmax,Psgmin);
  
    time(iter,1)=iteration-1;
    theta(:,iter)=theta_0;
    w(:,iter)=omega_0;
    Edp(:,iter)= E_dp0;
    Eqp(:,iter)= E_qp0;
    Efd(:,iter)= E_fd0;
    Pm(:,iter)= Pm_0;
    iter=iter+1;
end
iter=iter-1; %% To delete the extra time step added at the end
%     plot(time,rad2deg(theta(2,:)),'linewidth',1.5)

%% plot graphs
 
figure(1)
plot(time,theta(2:4,:)); grid on;
ylabel('\theta');
xlabel('Time/s');
legend('Gen2','Gen3','Gen4' )

figure(2)
plot(time,w(2:4,:)); grid on;
ylabel('\omega');
xlabel('Time/s');
legend('Gen2','Gen3','Gen4' )

figure(3)
plot(time,Edp(2:4,:)); grid on;
ylabel('Edp');
xlabel('Time/s');
legend('Gen2','Gen3','Gen4' )

figure(4)
plot(time,Eqp(2:4,:)); grid on;
ylabel('Eqp');
xlabel('Time/s');
legend('Gen2','Gen3','Gen4' )

figure(5)
plot(time,Efd(2:4,:)); grid on;
ylabel('Efd');
xlabel('Time/s');
legend('Gen2','Gen3','Gen4' )

figure(6)
plot(time,Pm(2:4,:)); grid on;
ylabel('Pm');
xlabel('Time/s');
legend('Gen2','Gen3','Gen4' )






