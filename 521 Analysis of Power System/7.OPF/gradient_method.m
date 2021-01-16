%% gradient_method
function [gradient_L_u,P12]=gradient_method(P1,P2,P12_max,J,V_mag,V_Delta,Theta,Y_mag,PQ,nPQ)

gradient_f_u = 0.0096*P2*100+6.4;
gradient_g_u=[1;zeros(21,1)];
gradient_g_x=-J;

gradient_P1_x=[];


for i=2:14
    gradient_P1_x=[gradient_P1_x;V_mag(1)*V_mag(i)*Y_mag(1,i)*sin(V_Delta(1)-V_Delta(i)-Theta(1,i))];
end
for j=1:nPQ
    gradient_P1_x=[gradient_P1_x;V_mag(1)*Y_mag(1,PQ(j))*cos(V_Delta(1)-V_Delta(PQ(j))-Theta(1,PQ(j)))];
end
P12 =[V_mag(1)*V_mag(2)*Y_mag(1,2)*cos(V_Delta(1)-V_Delta(2)-Theta(1,2))-V_mag(1)^2*Y_mag(1,2)*cos(Theta(1,2))]*100;

if abs(P12)>P12_max
    gradient_P12_x=[V_mag(1)*V_mag(2)*Y_mag(1,2)*sin(V_Delta(1)-V_Delta(2)-Theta(1,2));zeros(21,1)];
else
    gradient_P12_x=0;
end
    
gradient_f_x=(0.008*P1+8)*gradient_P1_x+2*(P12-P12_max)*gradient_P12_x;

gradient_L_u=gradient_f_u-gradient_g_u'*inv(gradient_g_x')*gradient_f_x;
    
        

