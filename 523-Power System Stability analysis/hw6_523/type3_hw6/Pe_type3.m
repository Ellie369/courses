function [Pe]=Pe_type3(Y_gen_mag,Y_gen_angle,E_qp_mag,E_qp_angle,nPV,PV)

Pe=zeros(4,1);
for m=1:nPV
    i=PV(m);
    for j=1:4        % include slack bus
        Pe(i)= Pe(i)+Y_gen_mag(i,j)*E_qp_mag(i)*E_qp_mag(j)*cos(E_qp_angle(i)-E_qp_angle(j)-Y_gen_angle(i,j))        
    end   
end
end