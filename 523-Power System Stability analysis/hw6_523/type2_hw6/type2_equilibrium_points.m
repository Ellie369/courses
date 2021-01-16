%% x1---theta    x2--omega  
%%   x3---Eq'       x4---Ed'  
%%  x5---Id       x6---Iq   
%%  x7---Efd        x8---Pm  
%%  x9--- Vref       x10---Pc 
%%  x11--- pe     

function [f]=type2_equilibrium_points(x,genbus,ws,K_D,X_d,X_dp,X_q,X_qp,KA,P_gen_cal, Q_gen_cal, V_bus0,Delta0,R) 

% x(11)= P_gen_cal(genbus);
% x(10)= P_gen_cal(genbus);
% x(8)= P_gen_cal(genbus);

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
% f7_1=0;
% for j=1:4 % all generators, including slack bus
%     f7_1=f7_1+(Y_gen_mag(genbus,j)*E_p_mag(genbus)*E_p_mag(j)*cos(E_p_angle(genbus)-E_p_angle(j)-Y_gen_angle(genbus,j))); 
% end 
% f7=-x(11)+ f7_1; % Pe=...

% f8_1=0;
% for j=1:4    
%     f8_1=f8_1+(Y_gen_mag(genbus,j)*E_p_mag(j)*cos(x(1)-E_p_angle(j)-Y_gen_angle(genbus,j)));
% end
% f8= -x(6)+ f8_1; % I_q...
% 
% f9_1=0;
% for j=1:4    
%     f9_1=f9_1+(Y_gen_mag(genbus,j)*E_p_mag(j)*sin(x(1)-E_p_angle(j)-Y_gen_angle(genbus,j)));
% end
% f9= -x(5)+ f9_1; %I_d...
% 
% f10= -x(3)+E_p_mag(genbus)*cos(x(1)-E_p_angle(genbus)); %E_qp=E_p_mag*cos(theta-gama)
% f11 = -E_p_angle(genbus) +atan(x(3)/x(4))+x(1)- pi/2; %gama=arctan(E_qp/E_E_dp)+theta-pi/2;
% f12= -x(4)+E_p_mag(genbus)*sin(x(1)-E_p_angle(genbus)); % E_dp=E_p_mag*sin(theta-gama)

f=[f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11];

% % x1---theta    x2--omega  
% %%   x3---Eq'       x4---Ed'  
% %%  x5---Id       x6---Iq   
% %%  x7---Efd        x8---Pm  
% %%  x9--- Vref       x10---Pc 
% %%  x11--- pe 
% 
% function [f]=type2_equilibrium_points(x,genbus,P_gen_cal,ws,K_D,X_d,X_dp,X_q,X_qp,KA,V_mag,V_Delta,R) 
% f1=(x(2)-1)*ws; %theta_dot=(omega-1)*omega_s
% f2= x(8)-x(11)-K_D(genbus)*(x(2)-1); % 2H*omega_dot=Pm-Pe-Kd(omega-1)
% f3= (-x(3)-(X_d(genbus)-X_dp(genbus))*x(5)+x(7)); %Tdo_p*Eqp_dot=-E_qp-(X_d-X_dp)*I_d+Efd
% f4=(-x(4)+(X_q(genbus)-X_qp(genbus))*x(6)); %Tq0_p*Edp_dot=-E_dp+(X_q-X_qp)*I_q
% f5=(-x(7)+(x(9)-V_mag(genbus))*KA); %exciter: V_Ai_dot*T_Ai=(-V_Ai+K_Ai(V_refi-Vi))
% f6=(-x(8)+x(10)+(1/R)*(1-x(2))); %governer: Psgi_dot*Tsgi=(-Psgi+Ksgi(Pci+(1/Ri)*(ws-wi))),usually Ksgi=1 
% 
% % E_p_mag(genbus)= sqrt(x(3)^2+x(4)^2);
% % E_p_angle(genbus)= atan(x(3)/x(4));
% 
% 
% % f7_1=0;
% % for j=1:4 % all generators, including slack bus
% %     f7_1=f7_1+(Y_gen_mag(genbus,j)*E_p_mag(genbus)*E_p_mag(j)*cos(E_p_angle(genbus)-E_p_angle(j)-Y_gen_angle(genbus,j))); 
% % end 
% % f7=-x(11)+ f7_1; % Pe=...
% f7= -x(11)+P_gen_cal;
% % f8_1=0;
% % for j=1:4    
% %     f8_1=f8_1+(Y_gen_mag(genbus,j)*E_p_mag(j)*cos(x(1)-E_p_angle(j)-Y_gen_angle(genbus,j)));
% % end
% % f8= -x(6)+ f8_1; % I_q...
% f8=-x(6)+( -1/X_qp(genbus))* (x(4)-V_mag(genbus)*sin(x(1)-V_Delta(genbus)));
% 
% 
% % f9_1=0;
% % for j=1:4    
% %     f9_1=f9_1+(Y_gen_mag(genbus,j)*E_p_mag(j)*sin(x(1)-E_p_angle(j)-Y_gen_angle(genbus,j)));
% % end
% % f9= -x(5)+ f9_1; %I_d...
% f9=-x(5)+( 1/X_dp(genbus))* (x(3)-V_mag(genbus)*cos(x(1)-V_Delta(genbus)));
% 
% % f10= -x(3)+E_p_mag(genbus)*cos(x(1)-E_p_angle(genbus)); %E_qp=E_p_mag*cos(theta-gama)
% % f11 = -E_p_angle(genbus) +atan(x(3)/x(4))+x(1)- pi/2; %gama=arctan(E_qp/E_E_dp)+theta-pi/2;
% % f12= -x(4)+E_p_mag(genbus)*sin(x(1)-E_p_angle(genbus)); % E_dp=E_p_mag*sin(theta-gama)
% f10= -x(8)+ P_gen_cal;
% f11= -x(10)+ P_gen_cal;
% 
% f=[f1;f2;f3;f4;f5;f6;f7;f8;f9;f10;f11];
% end






