%% homework 6
%generator parameters, changing rating from 900MVA to 100MVA 
%p---prime; pp---double prime
clc
clear all
 format shortEng
 format compact
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
R=0.05/9;
%% read grid data from file
file_name='b_kundur_system.txt';
[S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data]=read_data(file_name);
[PQ,nPQ,PV,nPV,Y_ang,Y_mag,V_mag,V_Delta,P_gen_cal,Q_gen_cal, V_result, P_load, Q_load]= NR_power_flow_type1(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data);

%% solve equilibrium points 
V_bus0=V_mag;
Delta0=V_Delta;
for genbus=1:4 %generator buses, including slack bus
    x0=[0;1;1;0.5;1;1;1.7;7;1;7;7]; %
    options=optimoptions('fsolve','algorithm','levenberg-marquardt','display','off');
    [solution,~,exitflag]  = fsolve(@(x) type1_equilibrium_points(x,genbus,ws,K_D,X_d,X_dp,X_q,X_qp,KA,P_gen_cal, Q_gen_cal, V_bus0,Delta0,R),x0,options);
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
%% consider load as ZIP model
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


V_bus=sym('V_bus', [1,11]);
% V_bus_mag= sym('V_bus_mag', [1,11]);
Delta=sym('Delta',[1,11]);
theta=sym('theta',[1,4]); %rotor angle
omega=sym('omega',[1,4]);
Eqp=sym('Eqp',[1,4]);
Edp=sym('Edp',[1,4]);
Efd=sym('Efd',[1,4]);
Psg=sym('Psg',[1,4]);

Pg= zeros(1,11)*omega(1);
Qg= zeros(1,11)*omega(1);

V_bus_mag=sqrt(real(V_bus).^2+imag(V_bus).^2);
for i=1:4
Vd(i)=V_bus_mag(i)*sin(theta(i)-Delta(i));
Vq(i)=V_bus_mag(i)*cos(theta(i)-Delta(i));

Id(i)=(Eqp(i)-V_bus_mag(i)*cos(theta(i)-Delta(i)))/X_dp(i);
Iq(i)=(-1/X_qp(i))*(Edp(i)-V_bus_mag(i)*sin(theta(i)-Delta(i)));
Pg(i)=Vd(i)*Id(i)+Vq(i)*Iq(i);
Qg(i)=Vq(i)*Id(i)-Vd(i)*Iq(i);
end

for i=2:4

    f(1+6*(i-2))=(omega(i)-1)*ws;
    f(2+6*(i-2))=(Psg(i)-Pg(i)-K_D(i)*(omega(i)-1))/(2*H(i));
    f(3+6*(i-2))=(-Eqp(i)-(X_d(i)-X_dp(i))*Id(i)+Efd(i))/T_d0p(i);
    f(4+6*(i-2))=(-Edp(i)+(X_q(i)-X_qp(i))*Iq(i))/T_q0p(i);
    f(5+6*(i-2))=(-Efd(i)+KA*(V_ref(i)-V_bus_mag(i)))/TA;
    f(6+6*(i-2))=(-Psg(i)+Ksg*(Pc(i)+(1-omega(i))/R))/Tsg;

end


x0_A=[theta_0(2:4)',omega_0(2:4)',E_qp0(2:4)',E_dp0(2:4)',E_fd0(2:4)',Pm_0(2:4)',Delta0(2:4)',V_bus0(2:4)'];
 A=jacobian(f,[theta(2),omega(2),Eqp(2),Edp(2),Efd(2),Psg(2),...
    theta(3),omega(3),Eqp(3),Edp(3),Efd(3),Psg(3),...
    theta(4),omega(4),Eqp(4),Edp(4),Efd(4),Psg(4)]);
A=double(subs(A,[theta(2:4),omega(2:4),Eqp(2:4),Edp(2:4),Efd(2:4),Psg(2:4),Delta(2:4),V_bus(2:4)],x0_A));
x0_B=[theta_0(2:4)',omega_0(2:4)',E_qp0(2:4)',E_dp0(2:4)',E_fd0(2:4)',Pm_0(2:4)',Delta0(2:11)',V_bus0(2:11)'];
B= jacobian(f,[V_bus(2:11),Delta(2:11)]);
B= double(subs(B,[theta(2:4),omega(2:4),Eqp(2:4),Edp(2:4),Efd(2:4),Psg(2:4),Delta(2:11),V_bus(2:11)],x0_B));

P_load=PL0'+(M0').*V_bus_mag + (G0').*(V_bus_mag.^2);
Q_load=QL0'+(H0').*V_bus_mag + (B0').*(V_bus_mag.^2);
P_temp= zeros(1,11)*omega(1);
Q_temp= zeros(1,11)*omega(1);
for i=2:11
    for j=1:11
    P_temp(i)=P_temp(i)+V_bus_mag(i)*V_bus_mag(j)*Y_mag(i,j)*cos(Delta(i)-Delta(j)-Y_ang(i,j));
    Q_temp(i)=Q_temp(i)+V_bus_mag(i)*V_bus_mag(j)*Y_mag(i,j)*sin(Delta(i)-Delta(j)-Y_ang(i,j));
    end    
    g(i-1)=Pg(i)-P_load(i)-P_temp(i);
    g(i+9)=Qg(i)-Q_load(i)-Q_temp(i);   
end
x0_C=[theta_0(2:4)',omega_0(2:4)',E_qp0(2:4)',E_dp0(2:4)',E_fd0(2:4)',Pm_0(2:4)',Delta0(1:11)',V_bus0(1:11)'];
C= jacobian(g,[theta(2),omega(2),Eqp(2),Edp(2),Efd(2),Psg(2),...
    theta(3),omega(3),Eqp(3),Edp(3),Efd(3),Psg(3),...
    theta(4),omega(4),Eqp(4),Edp(4),Efd(4),Psg(4)]);
C= double(subs(C,[theta(2:4),omega(2:4),Eqp(2:4),Edp(2:4),Efd(2:4),Psg(2:4),Delta(1:11),V_bus(1:11)],x0_C))
 
x0_D=[theta_0(2:4)',omega_0(2:4)',E_qp0(2:4)',E_dp0(2:4)',E_fd0(2:4)',Pm_0(2:4)',Delta0(1:11)',V_bus0(1:11)'];
D=jacobian(g,[V_bus(2:11),Delta(2:11)]);
D= double(subs(D,[theta(2:4),omega(2:4),Eqp(2:4),Edp(2:4),Efd(2:4),Psg(2:4),Delta(1:11),V_bus(1:11)],x0_D));

J_type1=A-B*inv(D)*C;
[right_EV,Eigen]=eig(J_type1);     %% Right eigen vector and eigenvalues
left_EV=inv(right_EV);              %% Left eigen vector

EG=eig(J_type1)
plot(EG,'o')
axis([-100 5 -8 8])
xlabel('Real')
ylabel('Imaginary')
title('Eigenvalue for Type 1 model')

%% participation factor matrix
for i=1:length(J_type1)
    for k=1:length(J_type1)
        Participation_matrix(k,i)=right_EV(k,i)*left_EV(i,k);
    end        
end

Max_Participation_matrix=max(abs(Participation_matrix)); % Normalaizing the Participation matrix by dividing by the maximum valueof each column

for e=1:length(J_type1)
    Participation_matrix(:,e)=abs(Participation_matrix(:,e))/Max_Participation_matrix(e);
end
%% Calculating the Frequency
frequency = abs(imag(EG))./(2*pi)

%% Calculating the damping ratio
damping_ratio = -real(EG)./abs(EG)



