%% homework 6
%generator parameters, changing rating from 900MVA to 100MVA 
%p---prime; pp---double prime
clc
clear all
format short 
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
[PQ,nPQ,PV,nPV,Y_mat,V_mag,V_Delta,P_gen_cal,Q_gen_cal, V_result]= NR_power_flow(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data)

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
E_p(i)=complex(V_mag(i)*cos(V_Delta(i)),V_mag(i).*sin(V_Delta(i)))+I(i).*(R_a(i)+1i.*X_p(i))
end
E_p_mag=abs(E_p);
E_p_angle=angle(E_p);
w_0(PV)=1; %initialize omega
[Pe_0]=Pe_type3(Y_gen_mag,Y_gen_angle,E_p_mag,E_p_angle,nPV,PV);
Pm_0=Pe_0;

%% calculate Jacobian for model analysis
J_type3=[]
for e=1:nPV
     for k=1:nPV
         i=PV(e);
         j=PV(k);
        [A]=Jacobian_type3(i,j,ws,Y_gen_mag,Y_gen_angle,E_p_mag,E_p_angle,H,K_D)
        row = no_of_states*(e-1)+1;
        col = no_of_states*(k-1)+1;
        J_type3((row:(row+no_of_states-1)),(col:(col+no_of_states-1))) = A;  %% Store in the Jacobian structure
    end      
end

[right_EV,Eigen]=eig(J_type3);     %% Right eigen vector and eigenvalues
left_EV=inv(right_EV);              %% Left eigen vector

EG=eig(J_type3)
plot(EG,'r*')
axis([-1 1 -8 8])
xlabel('Real')
ylabel('Imaginary')
title('Eigenvalue for Type 3 model')
%% participation factor matrix
for i=1:length(J_type3)
    for k=1:length(J_type3)
        Participation_matrix(k,i)=right_EV(k,i)*left_EV(i,k);
    end        
end
% Normalaizing the Participation matrix by dividing by the maximum valueof each column

Max_Participation_matrix=max(abs(Participation_matrix));

for e=1:length(J_type3)
    Participation_matrix(:,e)=abs(Participation_matrix(:,e))/Max_Participation_matrix(e);
end

%% Calculating the Frequency
frequency = abs(imag(EG))./(2*pi)

%% Calculating the damping ratio
damping_ratio = -real(EG)./abs(EG)



