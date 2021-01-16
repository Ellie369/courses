%% x1---theta    x2--omega  
%%   x3---Eq'       x4---Ed'  
%%  x5---Id       x6---Iq   
%%  x7---Efd        x8---Pm  
%%  x9--- Vref       x10---Pc 
%%  x11--- pe     

function [f]=type1_equilibrium_points(x,genbus,ws,K_D,X_d,X_dp,X_q,X_qp,KA,P_gen_cal, Q_gen_cal, V_bus0,Delta0,R) 

f1=(x(2)-1)*ws; %theta_dot=(omega-1)*omega_s
f2= x(8)-x(11)-K_D(genbus)*(x(2)-1); % 2H*omega_dot=Pm-Pe-Kd(omega-1)
f3= (-x(3)-(X_d(genbus)-X_dp(genbus))*x(5)+x(7)); %Tdo_p*Eqp_dot=-E_qp-(X_d-X_dp)*I_d+Efd
f4=(-x(4)+(X_q(genbus)-X_qp(genbus))*x(6)); %Tq0_p*Edp_dot=-E_dp+(X_q-X_qp)*I_q
f5=(-x(7)+(x(9)-V_bus0(genbus))*KA); %exciter: V_Ai_dot*T_Ai=(-V_Ai+K_Ai(V_refi-Vi))
f6=(-x(8)+x(10)+(1/R)*(1-x(2))); %governer: Psgi_dot*Tsgi=(-Psgi+Ksgi(Pci+(1/Ri)*(ws-wi))),usually Ksgi=1 


Vd(genbus)=V_bus0(genbus)*sin(x(1)-Delta0(genbus)); %Vdi=Vi*sin(theta-delta)
Vq(genbus)=V_bus0(genbus)*cos(x(1)-Delta0(genbus)); %Vqi=Vi*cos(theta-delta)

f7=-P_gen_cal(genbus)+Vd(genbus)*x(5)+Vq(genbus)*x(6); % Pgi=Vdi*Idi+Vqi*Iqi
f8=-Q_gen_cal(genbus)+Vq(genbus)*x(5)-Vd(genbus)*x(6); % Qgi= Vqi*Idi-Vdi*Iqi

f9=-x(3)+Vq(genbus)+X_dp(genbus)*x(5); % Eq'=Vqi+Xd'Idi+RsiIqi
f10=-x(4)+Vd(genbus)-X_qp(genbus)*x(6); % Ed'=Vdi-Xq'Iqi+RsiIdi

f11=-x(11)+P_gen_cal(genbus);

f=[f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11];
end






