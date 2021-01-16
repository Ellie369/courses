%% homework 6
%generator parameters, changing rating from 900MVA to 100MVA 
%p---prime; pp---double prime
clc
clear all
format shortEng
format compact
%% Type2 model fault
 
fault_time =2;    %%t=t+1   %% This is the time when the fault occurs
clearing_time = 3/60;        %% This is the time of clearance after the fault has occurred
endtime = 51;    %% t=t+1
h = 1/1000;   %% Step size for simulation


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
% for type2
X_dp=X_p;
X_qp=X_p;
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
Psgmax=1*9;
R=0.0056; %convert to Base =100MVA
%% read grid data from file
file_name='b_kundur_system.txt';
[S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data]=read_data(file_name);
[PQ,nPQ,PV,nPV,Y_mat,V_mag,V_Delta,P_gen_cal,Q_gen_cal, V_result]= NR_power_flow(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data);

%% revise y matrix, and form y_gen
file_name_2='b_kundur_system_extended.txt';
[S_Base_2,No_of_Buses_2,No_of_Lines_2,Bus_data_2,Line_data_2]=read_data(file_name_2);
[Y_mat_ex,Theta_ex,Y_mag_ex,B_ex,G_ex]=y_bus(Bus_data_2,Line_data_2,No_of_Buses_2,No_of_Lines_2);
[Y_gen,Y_gen_mag,Y_gen_angle]=Ygen(Y_mat_ex,4,10);

%% solving for equilibrium point at t=0
S=complex(P_gen_cal, Q_gen_cal);
I=conj(S./(V_mag.*cos(V_Delta)+1i.*V_mag.*sin(V_Delta)));
I_mag=abs(I);
I_angle=angle(I);
E_p=zeros(4,1);
for i=1:4
E_p(i)=complex(V_mag(i)*cos(V_Delta(i)),V_mag(i).*sin(V_Delta(i)))+I(i).*(R_a(i)+1i.*X_p(i));
end
E_p_mag=abs(E_p); %E'
E_p_angle=angle(E_p); %gama

for genbus=1:4 %generator buses, including slack bus
     x0=[0;1;1.06;0;1;1;1;1;1;1;1]; %
    [solution,~,exitflag] = fsolve(@(x) type2_equilibrium_points(x,genbus,ws,K_D,X,X_p,X,X_p,KA,V_mag,Y_gen_mag,Y_gen_angle,E_p_mag,E_p_angle,R),x0,optimset('algorithm','levenberg-marquardt','display','off'));
    theta_0(genbus,1)=solution(1);
    omega_0(genbus,1)=solution(2);
    E_qp0(genbus,1)=solution(3);
    E_dp0(genbus,1)=solution(4);
    E_fd0(genbus,1)=solution(7);
    Pm_0(genbus,1)=solution(8);
    V_ref(genbus,1)=solution(9);
    Pc_0(genbus,1)=solution(10);
    Pe_0(genbus,1)=solution(11);
    Id_0(genbus,1)=solution(5);
    Iq_0(genbus,1)=solution(6);
     
end

%% Transient stability analysis
%% prefault
iter=1;
for iteration=1:h:fault_time
   [theta_0,omega_0,E_dp0, E_qp0,E_fd0, Pm_0,V_mag]=type2_euler(Pm_0, Pc_0,PV, nPV,H,Y_gen,Y_gen_mag,Y_gen_angle,V_mag,theta_0,omega_0,ws,K_D,h,...
        E_qp0, E_dp0, E_fd0, V_ref, X_d, X_q, X_dp, X_qp,X_p, T_d0p,T_q0p,KA,TA,Tsg,Ksg,R, VRmax,VRmin,Efdmax,Efdmin,Psgmax,Psgmin);
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
%% onfault
% revise y matrix, and form y_gen
file_name_3='onfault_b_kundur_system_extended.txt';
[S_Base_3,No_of_Buses_3,No_of_Lines_3,Bus_data_3,Line_data_3]=read_data(file_name_3);
[Y_mat_ex_f,Theta_ex_f,Y_mag_ex_f,B_ex_f,G_ex_f]=y_bus(Bus_data_3,Line_data_3,No_of_Buses_3,No_of_Lines_3);
[Y_gen_f,Y_gen_mag_f,Y_gen_angle_f]=Ygen(Y_mat_ex_f,4,11);

for iteration=fault_time:h:(fault_time+clearing_time)
     [theta_0,omega_0,E_dp0, E_qp0,E_fd0, Pm_0,V_mag]=type2_euler(Pm_0, Pc_0,PV, nPV,H,Y_gen_f,Y_gen_mag_f,Y_gen_angle_f,V_mag,theta_0,omega_0,ws,K_D,h,...
        E_qp0, E_dp0, E_fd0, V_ref, X_d, X_q, X_dp, X_qp,X_p, T_d0p,T_q0p,KA,TA,Tsg,Ksg,R,VRmax,VRmin,Efdmax,Efdmin,Psgmax,Psgmin);
    time(iter,1)=iteration-1;
    theta(:,iter)=theta_0;
    w(:,iter)=omega_0;
    Edp(:,iter)= E_dp0;
    Eqp(:,iter)= E_qp0;
    Efd(:,iter)= E_fd0;
    Pm(:,iter)= Pm_0;
iter = iter+1;  
end
iter=iter-1; %% To delete the extra time step added at the end

%% post fault analysis
% revise y matrix, and form y_gen
file_name_4='postfault_b_kundur_system_extended.txt';
[S_Base_4,No_of_Buses_4,No_of_Lines_4,Bus_data_4,Line_data_4]=read_data(file_name_4);
[Y_mat_ex_pf,Theta_ex_pf,Y_mag_ex_pf,B_ex_pf,G_ex_pf]=y_bus(Bus_data_4,Line_data_4,No_of_Buses_4,No_of_Lines_4);
[Y_gen_pf,Y_gen_mag_pf,Y_gen_angle_pf]=Ygen(Y_mat_ex_pf,4,10);

for iteration=(fault_time+clearing_time):h:endtime
   [theta_0,omega_0,E_dp0, E_qp0,E_fd0, Pm_0,V_mag]=type2_euler(Pm_0, Pc_0,PV, nPV,H,Y_gen_pf,Y_gen_mag_pf,Y_gen_angle_pf,V_mag,theta_0,omega_0,ws,K_D,h,...
        E_qp0, E_dp0, E_fd0, V_ref, X_d, X_q, X_dp, X_qp,X_p, T_d0p,T_q0p,KA,TA,Tsg,Ksg,R,VRmax,VRmin,Efdmax,Efdmin,Psgmax,Psgmin);
    time(iter,1)=iteration-1;
    theta(:,iter)=theta_0;
    w(:,iter)=omega_0;
    Edp(:,iter)= E_dp0;
    Eqp(:,iter)= E_qp0;
    Efd(:,iter)= E_fd0;
    Pm(:,iter)= Pm_0;
iter = iter+1;  
end
iter=iter-1; %% To delete the extra time step added at the end

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







