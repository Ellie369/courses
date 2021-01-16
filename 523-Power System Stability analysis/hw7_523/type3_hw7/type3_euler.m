function [theta_0,w_0]=type3_euler(Pe_0, Pm_0, PV, nPV,H,Y_gen_mag,Y_gen_angle,E_qp_mag,theta_0,w_0,ws,K_D,h)
% calculate Pe
[Pe]=Pe_type3(Y_gen_mag,Y_gen_angle,E_qp_mag,theta_0,nPV,PV);
% calculate w_1
for m=1:nPV
     i=PV(m);
     w_dot(i,1)=1/(2*H(i))*(Pm_0(i)-Pe(i)-K_D(i)*(w_0(i)-1));
     w_0(i)=w_0(i)+h.*w_dot(i);
end


% calculate theta_0
for m=1:nPV
    i=PV(m);
    theta_dot(i,1)=(w_0(i)-1)*ws;    
    theta_0(i)=theta_0(i)+h.*theta_dot(i);
end





    


