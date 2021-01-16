function [A]=Jacobian_type3(i,j,ws,Y_gen_mag,Y_gen_angle,E_qp_mag,E_qp_angle,nPV,PV,H,K_D)
if i==j
    A11=0;
    A12=ws;
    A21=0;
    for m=1:4        
        if i~=m
        A21=A21+Y_gen_mag(i,m)*E_qp_mag(i)*E_qp_mag(m)*sin(E_qp_angle(i)-E_qp_angle(m)-Y_gen_angle(i,m))    
        end
    end 
    A21=1/(2*H(i))*A21;
    A22=-K_D(i)/(2*H(i));
end
if i~=j
    A11=0;
    A12=0;
    A21=-1/(2*H(i))*(Y_gen_mag(i,j)*E_qp_mag(i)*E_qp_mag(j)*sin(E_qp_angle(i)-E_qp_angle(j)-Y_gen_angle(i,j)))        
    A22=0;
end
A=[A11,A12;A21,A22];
end