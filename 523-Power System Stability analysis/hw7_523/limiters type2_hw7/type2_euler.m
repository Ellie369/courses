function  [theta_0,omega_0,E_dp0, E_qp0,E_fd0, Pm_0,V_bus_mag]=type2_euler(Pm_0, Pc_0,PV, nPV,H,Y_gen,Y_gen_mag,Y_gen_angle,V_bus_mag,theta_0,omega_0,ws,K_D,h,...
        E_qp0, E_dp0, E_fd0, V_ref, X_d, X_q, X_dp, X_qp,X_p, T_d0p,T_q0p,KA,TA,Tsg,Ksg,R,VRmax,VRmin,Efdmax,Efdmin,Psgmax,Psgmin)
%% calculate Pe
% [Pe_0]=Pe_type2(Y_gen_mag,Y_gen_angle,E_qp0,E_dp0, theta_0);
for i=1:4
    E_p(i)=sqrt((E_dp0(i))^2+(E_qp0(i))^2)*exp(1j*(atan(E_qp0(i)/E_dp0(i))+theta_0(i)-pi/2));
    E_p_mag(i)=sqrt(E_dp0(i)^2+E_qp0(i)^2);
    gama(i)=atan(E_qp0(i)/E_dp0(i))+theta_0(i)-pi/2;
end
% update Pe, Id, Iq, V_bus_mag
Pe_0(2:4)= 0;

for i=2:4
    for j=1:4
        Pe_0(i)=Pe_0(i)+Y_gen_mag(i,j)*E_p_mag(i)*E_p_mag(j)*cos(gama(i)-gama(j)-Y_gen_angle(i,j));
    end
end

%% calculate omega_0
for j=1:nPV
     i=PV(j);
     w_dot(i,1)=1/(2*H(i))*(Pm_0(i)-Pe_0(i)-K_D(i)*(omega_0(i)-1));
     omega_0(i)=omega_0(i)+h.*w_dot(i);
end
% calculate theta_0
for j=1:nPV
    i=PV(j);
    theta_dot(i,1)=(omega_0(i)-1)*ws;    
    theta_0(i)=theta_0(i)+h.*theta_dot(i);
end
% update Id, Iq
Id_0(2:4)= 0;
Iq_0(2:4)= 0;
for i=2:4
    for j=1:4
        Id_0(i)=Id_0(i)+Y_gen_mag(i,j)*E_p_mag(j)*sin(theta_0(i)-gama(j)-Y_gen_angle(i,j));
        Iq_0(i)=Iq_0(i)+Y_gen_mag(i,j)*E_p_mag(j)*cos(theta_0(i)-gama(j)-Y_gen_angle(i,j));
    end
end
% calculate Eq'
for j=1:nPV
    i=PV(j);
    Eqp_dot(i,1)=(-E_qp0(i)-(X_d(i)-X_dp(i))*Id_0(i)+E_fd0(i))/T_d0p(i);
    E_qp0(i)=E_qp0(i)+h.*Eqp_dot(i);
end
% calculate Ed'
for j=1:nPV
    i=PV(j);
    Edp_dot(i,1)= (-E_dp0(i)+(X_q(i)-X_qp(i))*Iq_0(i))/T_q0p(i);
    E_dp0(i)=E_dp0(i)+h.*Edp_dot(i);
end

%% update Pe,Pe(1)--slack bus keep the same
% [Pe]=Pe_type2(Y_gen_mag,Y_gen_angle,E_qp_mag,theta_0,nPV,PV);
% for i=1:4
%     E_p(i)=sqrt((E_dp0(i))^2+(E_qp0(i))^2)*exp(1j*(atan(E_qp0(i)/E_dp0(i))+theta_0(i)-pi/2));
%     E_p_mag(i)=sqrt(E_dp0(i)^2+E_qp0(i)^2);
%     gama(i)=atan(E_qp0(i)/E_dp0(i))+theta_0(i)-pi/2;
% end
% % update Pe, Id, Iq, V_bus_mag
% Pe_0(2:4)= 0;
% 
% for i=2:4
%     for j=1:4
%         Pe_0(i)=Pe_0(i)+Y_gen_mag(i,j)*E_p_mag(i)*E_p_mag(j)*cos(gama(i)-gama(j)-Y_gen_angle(i,j));
%     end
% end
I_gen=zeros(4,1);
for i=2:4
    for j=1:4
        I_gen(i)=I_gen(i)+Y_gen(i,j)*E_p(j);
    end 
    V_bus(i)=E_p(i)-I_gen(i)*(1j*X_p(i));
    V_bus_mag(i)=sqrt(real(V_bus(i))^2+imag(V_bus(i))^2);
end

% update Efd
for j=1:nPV
    i=PV(j);
    Efd_dot(i,1)= (-E_fd0(i)+KA*(V_ref(i)-V_bus_mag(i)))/TA;    
   
    % Exciter non-windup limiter
   
    if (E_fd0(i)>= VRmax)&&(Efd_dot(i,1)>0)
        Efd_dot(i,1)=0;
    elseif (E_fd0(i)<= VRmin)&&(Efd_dot(i,1)<0)
        Efd_dot(i,1)=0;
    end
    E_fd0(i)=E_fd0(i)+h.*Efd_dot(i); 
    % Exciter limiter (Windup limiter)
    
    if E_fd0(i)>Efdmax
        E_fd0(i)=Efdmax;
    elseif E_fd0(i)<Efdmin
        E_fd0(i)=Efdmin;
    end
   
end
%update Psg
for j=1:nPV
    i=PV(j);
    Pm_dot(i,1)= (-Pm_0(i)+Ksg*(Pc_0(i)+(1-omega_0(i))/R))/Tsg;
    % governer limiter
    if (Pm_0(i)<=Psgmin)&&(Pm_dot(i,1)<0)
        Pm_dot(i,1)=0;
    elseif (Pm_0(i)>=Psgmax)&&(Pm_dot(i,1)>0)
        Pm_dot(i,1)=0;
    end        
    Pm_0(i)=Pm_0(i)+h.*Pm_dot(i);
end
end

        
    


