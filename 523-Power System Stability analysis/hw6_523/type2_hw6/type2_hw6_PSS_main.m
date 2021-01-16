%% homework 6
%generator parameters, changing rating from 900MVA to 100MVA 
%p---prime; pp---double prime
clc
clear all
 format shortEng
 
%  format compact
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
A_sat(2:4,1)=0.015; 
B_sat(2:4,1)=9.6;  
psi_t1(2:4,1)=0.9;
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
no_of_states = 8; 
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
R=0.05/9;
%% parameters for pss
Kstab=5; 
T1=0.06;
T2=0.02;
Tw=10;

%% read grid data from file
file_name='b_kundur_system.txt';
[S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data]=read_data(file_name);
[PQ,nPQ,PV,nPV,Y_mat,V_mag,V_Delta,P_gen_cal,Q_gen_cal, V_result]= NR_power_flow(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data);

%% revise y matrix, and form y_gen
file_name_2='b_kundur_system_extended.txt';
[S_Base_2,No_of_Buses_2,No_of_Lines_2,Bus_data_2,Line_data_2]=read_data(file_name_2);
[Y_mat_ex,Theta_ex,Y_mag_ex,B_ex,G_ex]=y_bus(Bus_data_2,Line_data_2,No_of_Buses_2,No_of_Lines_2);
[Y_gen,Y_gen_mag,Y_gen_angle]=Ygen(Y_mat_ex,4,10);
% yl7=(9.67-1j)/V_mag(7)^2;
% yl9=(17.67-1j)/V_mag(9)^2;


% solving for equilibrium point at t=0
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
    solution  = fsolve(@(x) type2_pss_equilibrium_points(x,genbus,ws,K_D,X,X_p,X,X_p,KA,V_mag,Y_gen_mag,Y_gen_angle,E_p_mag,E_p_angle,R),x0,optimset('algorithm','levenberg-marquardt','display','off'));
    theta_0(genbus,1)=solution(1);
    omega_0(genbus,1)=solution(2);
    E_qp0(genbus,1)=solution(3);
    E_dp0(genbus,1)=solution(4);
    E_fd0(genbus,1)=solution(7);
    Pm_0(genbus,1)=solution(8);
    V_ref(genbus,1)=solution(9);
    Pc(genbus,1)=solution(10);
    Pe_0(genbus,1)=solution(11);
end
% initial values for PSS 
  Vw_0=zeros(1,4);
  Vpss_0=zeros(1,4);
V_bus=sym('V_bus', [1,4]);
V_bus_mag= sym('V_bus_mag', [1,4]);
theta=sym('theta',[1,4]); %rotor angle
omega=sym('omega',[1,4]);
Eqp=sym('Eqp',[1,4]);
Edp=sym('Edp',[1,4]);
Efd=sym('Efd',[1,4]);
Psg=sym('Psg',[1,4]);
Pe=zeros(4,1)*omega(1);
Id=zeros(4,1)*omega(1);
Iq=zeros(4,1)*omega(1);
I_gen=zeros(4,1)*omega(1);
Vw= sym('Vw',[1,4]);
Vpss = sym('Vpss',[1,4]);

E_p= sym('E_p',[1,4]);
E_p_mag= sym('E_p_mag',[1,4]);
gama= sym('gama',[1,4]);
for i=1:4
    E_p(i)=sqrt((Edp(i))^2+(Eqp(i))^2)*exp(1j*(atan(Eqp(i)/Edp(i))+theta(i)-pi/2));
    E_p_mag(i)=sqrt(Edp(i)^2+Eqp(i)^2);
    gama(i)=atan(Eqp(i)/Edp(i))+theta(i)-pi/2;
end

for i=1:4
    for j=1:4
        Pe(i)=Pe(i)+Y_gen_mag(i,j)*E_p_mag(i)*E_p_mag(j)*cos(gama(i)-gama(j)-Y_gen_angle(i,j));
        Id(i)=Id(i)+Y_gen_mag(i,j)*E_p_mag(j)*sin(theta(i)-gama(j)-Y_gen_angle(i,j));
        Iq(i)=Iq(i)+Y_gen_mag(i,j)*E_p_mag(j)*cos(theta(i)-gama(j)-Y_gen_angle(i,j));
        
    end
end

for i=1:4
    for m=1:4
    I_gen(i)=I_gen(i)+Y_gen(i,m)*E_p(m);
    end
% 
    V_bus(i)=E_p(i)-I_gen(i)*(1j*X_p(i));
   
    V_bus_mag(i)=sqrt(real(V_bus(i))^2+imag(V_bus(i))^2);
end
for i=2:4
    f(1+no_of_states*(i-2))=(omega(i)-1)*ws;
    f(2+no_of_states*(i-2))=(Psg(i)-Pe(i)-K_D(i)*(omega(i)-1))/(2*H(i));
    f(3+no_of_states*(i-2))=(-Eqp(i)-(X_d(i)-X_dp(i))*Id(i)+Efd(i))/T_d0p(i);
    f(4+no_of_states*(i-2))=(-Edp(i)+(X_q(i)-X_qp(i))*Iq(i))/T_q0p(i);
    f(5+no_of_states*(i-2))=(-Efd(i)+KA*(V_ref(i)-V_bus_mag(i)-Vpss(i)))/TA;
    f(6+no_of_states*(i-2))=(-Psg(i)+Ksg*(Pc(i)+(1-omega(i))/R))/Tsg;
    f(7+no_of_states*(i-2))=(1/Tw)*(-Vw(i)-f(2+no_of_states*(i-2))*Tw*Kstab);
    f(8+no_of_states*(i-2))=(1/T2)*(-Vpss(i)+Vw(i)+ f(7+no_of_states*(i-2))*T1);
   
end


x0_A=[theta_0(1:4)',omega_0(1:4)',E_qp0(1:4)',E_dp0(1:4)',E_fd0(1:4)',Pm_0(1:4)',Vw_0(1:4),Vpss_0(1:4)];
 J_type2=jacobian(f,[theta(2),omega(2),Eqp(2),Edp(2),Efd(2),Psg(2),Vw(2),Vpss(2),...
    theta(3),omega(3),Eqp(3),Edp(3),Efd(3),Psg(3),Vw(3),Vpss(3),...
    theta(4),omega(4),Eqp(4),Edp(4),Efd(4),Psg(4),Vw(4),Vpss(4)]);
J_type2=double(subs(J_type2,[theta(1:4),omega(1:4),Eqp(1:4),Edp(1:4),Efd(1:4),Psg(1:4),Vw(1:4),Vpss(1:4)],x0_A));


[right_EV,Eigen]=eig(J_type2);     %% Right eigen vector and eigenvalues
left_EV=inv(right_EV);              %% Left eigen vector

EG=eig(J_type2)

%% participation factor matrix
for i=1:length(J_type2)
    for k=1:length(J_type2)
        Participation_matrix(k,i)=right_EV(k,i)*left_EV(i,k);
    end        
end
% Normalaizing the Participation matrix by dividing by the maximum value of each column

Max_Participation_matrix=max(abs(Participation_matrix));

for e=1:length(J_type2)
    Participation_matrix(:,e)=Participation_matrix(:,e)/Max_Participation_matrix(e)
end

%% Calculating the Frequency
frequency = abs(imag(EG))./(2*pi)

%% Calculating the damping ratio
damping_ratio = -real(EG)./abs(EG)



