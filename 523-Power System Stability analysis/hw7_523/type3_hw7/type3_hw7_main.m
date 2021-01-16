%% homework 7
%generator parameters, changing rating from 900MVA to 100MVA 
%p---prime; pp---double prime
clc
clear all
format short 
%% Type3 model fault
 
fault_time = 11;    %%t=t+1   %% This is the time when the fault occurs
clearing_time = 3/60;        %% This is the time of clearance after the fault has occurred
endtime = 101;    %% t=t+1
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
A_sat(2:4,1)=0.015; % what's this
B_sat(2:4,1)=9.6;   % what's this
psi_t1(2:4,1)=0.9;
H(1)=6.5*900/100;
H(2)=6.5*900/100;
H(3)=6.175*900/100;
H(4)=6.175*900/100;
K_D(2:4,1)=2*900/100;

X=(X_d+X_q)/2;
X_p=(X_dp+X_qp)/2;
ws=2*pi*60;
no_of_states = 2; 

%% read grid data from file
file_name='b_kundur_system.txt';
[S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data]=read_data(file_name);
[PQ,nPQ,PV,nPV,Y_mat,V_mag,V_Delta,P_gen_cal,Q_gen_cal, V_result]= NR_power_flow(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data);

%% revise y matrix, and form y_gen
file_name_2='b_kundur_system_extended.txt';
[S_Base_2,No_of_Buses_2,No_of_Lines_2,Bus_data_2,Line_data_2]=read_data(file_name_2);
[Y_mat_ex,Theta_ex,Y_mag_ex,B_ex,G_ex]=y_bus(Bus_data_2,Line_data_2,No_of_Buses_2,No_of_Lines_2);

% yl7=(9.67-1j)/V_mag(7)^2;
% yl9=(17.67-1j)/V_mag(9)^2;
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
% Y_gen_mag=abs(Y_gen);
% Y_gen_angle=angle(Y_gen); %returns the phase angles in radians
E_p_mag=abs(E_p);
E_p_angle=angle(E_p);
theta_0=E_p_angle;
w_0(1:4,1)=1; %initialize omega
[Pe_0]=Pe_type3(Y_gen_mag,Y_gen_angle,E_p_mag,theta_0,nPV,PV);
Pm_0=Pe_0;

%% prefault
iter=1;
for iteration=1:h:fault_time
    [theta_0,w_0]=type3_euler(Pe_0, Pm_0, PV, nPV,H,Y_gen_mag,Y_gen_angle,E_p_mag,theta_0,w_0,ws,K_D,h);
    time(iter,1)=iteration-1;
    theta(:,iter)=theta_0;
    w(:,iter)=w_0;
    Pm(:,iter)=Pm_0;
    Vplot(:,iter)=V_mag(:,1);
    deltaplot(:,iter)=V_Delta(:,1);
    iter=iter+1;
end
iter=iter-1; %% To delete the extra time step added at the end

%% onfault
% revise y matrix, and form y_gen
file_name_3='onfault_b_kundur_system_extended.txt';
[S_Base_3,No_of_Buses_3,No_of_Lines_3,Bus_data_3,Line_data_3]=read_data(file_name_3);
[Y_mat_ex_f,Theta_ex_f,Y_mag_ex_f,B_ex_f,G_ex_f]=y_bus(Bus_data_3,Line_data_3,No_of_Buses_3,No_of_Lines_3);
[Y_gen_f,Y_gen_mag_f,Y_gen_angle_f]=Ygen(Y_mat_ex_f,4,11);

for iteration=fault_time:h:(fault_time+clearing_time)
    [theta_0,w_0]=type3_euler(Pe_0, Pm_0, PV, nPV,H,Y_gen_mag_f,...
        Y_gen_angle_f,E_p_mag,theta_0,w_0,ws,K_D,h);
    time(iter,1)=iteration-1;
    theta(:,iter)=theta_0; 
    w(:,iter)=w_0;
    Pm(:,iter)=Pm_0;
    Vplot(:,iter)=V_mag(:,1);
deltaplot(:,iter)=V_Delta(:,1);
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
    [theta_0,w_0]=type3_euler(Pe_0, Pm_0, PV, nPV,H,Y_gen_mag_pf,...
        Y_gen_angle_pf,E_p_mag,theta_0,w_0,ws,K_D,h);
    time(iter,1)=iteration-1;
    theta(:,iter)=theta_0; 
    w(:,iter)=w_0;
    Pm(:,iter)=Pm_0;
    Vplot(:,iter)=V_mag(:,1);
deltaplot(:,iter)=V_Delta(:,1);
iter = iter+1;  
end
iter=iter-1; %% To delete the extra time step added at the end

%% plot graphs

for e=1:nPV
    i=PV(e);    
figure(e)
subplot(2,1,1)
plot(time,rad2deg(theta(i,:)),'linewidth',1.5)
title((['Generator ' num2str(i) ' Type 3 model (clearing time = ' num2str(clearing_time*60) ' cycles)' ]),'FontSize',13);
ylabel('\theta','FontSize',12,'FontWeight','bold');
subplot(2,1,2)
plot(time,w(i,:),'linewidth',1.5)
ylabel('\omega','FontSize',14,'FontWeight','bold');

end






