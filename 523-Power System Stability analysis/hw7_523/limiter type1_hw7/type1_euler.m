function [theta_0,omega_0,E_dp0, E_qp0,E_fd0, Pm_0,P_gen,Q_gen,V_mag]=type1_euler(Pm_0, P_gen,Q_gen, Pc_0,H,V_mag,V_Delta, theta_0,omega_0,ws,K_D,h,...
        E_qp0, E_dp0, E_fd0, V_ref, X_d, X_q, X_dp, X_qp, T_d0p,T_q0p,KA,TA,Tsg,Ksg,R,VRmax,VRmin,Efdmax,Efdmin,Psgmax,Psgmin)
 %% calculate Pg, Id and Iq,Vd and Vq
for i=2:4
    Id(i) = (E_qp0(i)- V_mag(i)*cos(theta_0(i) -V_Delta(i)))/X_dp(i);
    Iq(i) = - (E_dp0(i) - V_mag(i)*sin(theta_0(i) - V_Delta(i)))/X_qp(i);  
    Vd(i)= V_mag(i)*sin(theta_0(i)-V_Delta(i));
    Vq(i)= V_mag(i)*cos(theta_0(i)-V_Delta(i));
end

for i=2:4   
    P_gen(i) = Vd(i)*Id(i)+Vq(i)* Iq(i);     
end
%% calculate omega_0

for i=2:4
     w_dot(i,1)=1/(2*H(i))*(Pm_0(i)-P_gen(i)-K_D(i)*(omega_0(i)-1));  
     omega_0(i)=omega_0(i)+h.*w_dot(i);
end
% calculate theta
for i=2:4
    theta_dot(i,1)=(omega_0(i)-1)*ws;    
    theta_0(i)=theta_0(i)+h.*theta_dot(i);
end

% calculate Eq'
for i=2:4
    Id(i) = (E_qp0(i)- V_mag(i)*cos(theta_0(i) -V_Delta(i)))/X_dp(i);
    Iq(i) = - (E_dp0(i) - V_mag(i)*sin(theta_0(i) - V_Delta(i)))/X_qp(i);  
end
for i=2:4
    Eqp_dot(i,1)=(-E_qp0(i)-(X_d(i)-X_dp(i))*Id(i)+E_fd0(i))/T_d0p(i);
    E_qp0(i)=E_qp0(i)+h.*Eqp_dot(i);
end
% calculate Ed'
for i=2:4
    Edp_dot(i,1)= (-E_dp0(i)+(X_q(i)-X_qp(i))*Iq(i))/T_q0p(i);
    E_dp0(i)=E_dp0(i)+h.*Edp_dot(i);
end
% update Efd
for i=2:4
   
    Vd(i)= V_mag(i)*sin(theta_0(i)-V_Delta(i));
    Vq(i)= V_mag(i)*cos(theta_0(i)-V_Delta(i));
    V_mag(i)=sqrt(Vd(i)^2+Vq(i)^2);
    P_gen(i) = Vd(i)*Id(i)+Vq(i)* Iq(i);
    Q_gen(i) = Vq(i)*Id(i)-Vd(i)* Iq(i);
   
end
for i=2:4
    Efd_dot(i,1)= (-E_fd0(i)+KA*(V_ref(i)-V_mag(i)))/TA;    
%      Exciter non-windup limiter
   
    if (E_fd0(i)>= VRmax)&&(Efd_dot(i,1)>0)
        Efd_dot(i,1)=0;
    elseif (E_fd0(i)<= VRmin)&&(Efd_dot(i,1)<0)
        Efd_dot(i,1)=0;
    end
    E_fd0(i)=E_fd0(i)+h.*Efd_dot(i); 
%     Exciter limiter (Windup limiter)
    
    if E_fd0(i)>Efdmax
        E_fd0(i)=Efdmax;
    elseif E_fd0(i)<Efdmin
        E_fd0(i)=Efdmin;
    end
   
    
end
%update Psg
for i=2:4
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

        
    


