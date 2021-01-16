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
psi_t1(2:4,1)=0.9;
H(1)=6.5*900/100;
H(2)=6.5*900/100;
H(3)=6.175*900/100;
H(4)=6.175*900/100;
K_D(2:4,1)=2*900/100;
ws=2*pi*60;
no_of_states = 6; 

% for type2 and type 3
X_p=(X_dp+X_qp)/2;
X_dp=X_p;
X_qp=X_p;

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
R=0.05/9
%% read grid data from file
file_name='b_kundur_system.txt';
[S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data]=read_data(file_name);
[PQ,nPQ,PV,nPV,Y_mat,V_mag,V_Delta,P_gen_cal,Q_gen_cal, V_result]= NR_power_flow(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data);

%% revise y matrix, and form y_gen
file_name_2='b_kundur_system_extended.txt';
[S_Base_2,No_of_Buses_2,No_of_Lines_2,Bus_data_2,Line_data_2]=read_data(file_name_2);
[Y_mat_ex,Theta_ex,Y_mag_ex,B_ex,G_ex]=y_bus(Bus_data_2,Line_data_2,No_of_Buses_2,No_of_Lines_2);
[Y_gen,Y_gen_mag,Y_gen_angle]=Ygen(Y_mat_ex,4,10);

V_bus0=V_mag;
Delta0=V_Delta;
for genbus=1:4 %generator buses, including slack bus
    x0=[0;1;1;0.5;1;1;1.7;7;1;7;7];  %
    [solution,~,exitflag]  = fsolve(@(x) type2_equilibrium_points(x,genbus,ws,K_D,X_d,X_dp,X_q,X_qp,KA,P_gen_cal, Q_gen_cal, V_bus0,Delta0,R) ,x0,optimset('algorithm','levenberg-marquardt','display','off'));
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

%% linearization

V_bus=sym('V_bus', [1,4]);
V_bus_mag= sym('V_bus_mag', [1,4]);
% Delta=sym('Delta',[1,11]);
theta=sym('theta',[1,4]); %rotor angle
omega=sym('omega',[1,4]);
Eqp=sym('Eqp',[1,4]);
Edp=sym('Edp',[1,4]);
Efd=sym('Efd',[1,4]);
Psg=sym('Psg',[1,4]);
Pe=zeros(4,1)*omega(1);
Id=zeros(4,1)*omega(1);
Iq=zeros(4,1)*omega(1);

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

I_gen=zeros(4,1)*omega(1);
for i=1:4
for m=1:4
    I_gen(i)=I_gen(i)+Y_gen(i,m)*E_p(m);
end 
    V_bus(i)=E_p(i)-I_gen(i)*(1j*X_p(i));
   
    V_bus_mag(i)=sqrt(real(V_bus(i))^2+imag(V_bus(i))^2);
end

for i=2:4

    f(1+6*(i-2))=(omega(i)-1)*ws;
    f(2+6*(i-2))=(Psg(i)-Pe(i)-K_D(i)*(omega(i)-1))/(2*H(i));
    f(3+6*(i-2))=(-Eqp(i)-(X_d(i)-X_dp(i))*Id(i)+Efd(i))/T_d0p(i);
    f(4+6*(i-2))=(-Edp(i)+(X_q(i)-X_qp(i))*Iq(i))/T_q0p(i);
    f(5+6*(i-2))=(-Efd(i)+KA*(V_ref(i)-V_bus_mag(i)))/TA;
    f(6+6*(i-2))=(-Psg(i)+Ksg*(Pc(i)+(1-omega(i))/R))/Tsg;

end


x0_A=[theta_0(1:4)',omega_0(1:4)',E_qp0(1:4)',E_dp0(1:4)',E_fd0(1:4)',Pm_0(1:4)'];
 J_type2=jacobian(f,[theta(2),omega(2),Eqp(2),Edp(2),Efd(2),Psg(2),...
    theta(3),omega(3),Eqp(3),Edp(3),Efd(3),Psg(3),...
    theta(4),omega(4),Eqp(4),Edp(4),Efd(4),Psg(4)]);
J_type2=double(subs(J_type2,[theta(1:4),omega(1:4),Eqp(1:4),Edp(1:4),Efd(1:4),Psg(1:4)],x0_A));


[right_EV,Eigen]=eig(J_type2);     %% Right eigen vector and eigenvalues
left_EV=inv(right_EV);              %% Left eigen vector

EG=eig(J_type2)
plot(EG,'b*')
axis([-100 5 -5 8])
xlabel('Real')
ylabel('Imaginary')
title('Eigenvalue for Type 2 model')


%% participation factor matrix
for i=1:length(J_type2)
    for k=1:length(J_type2)
        Participation_matrix(k,i)=right_EV(k,i)*left_EV(i,k);
    end        
end
% Normalaizing the Participation matrix by dividing by the maximum valueof each column

Max_Participation_matrix=max(abs(Participation_matrix));

for e=1:length(J_type2)
    Participation_matrix(:,e)=abs(Participation_matrix(:,e))/Max_Participation_matrix(e)
end
writematrix(Participation_matrix,'Participation2.csv')
%% Calculating the Frequency
frequency = abs(imag(EG))./(2*pi)

%% Calculating the damping ratio
damping_ratio = -real(EG)./abs(EG)



