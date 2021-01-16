function [Y_gen, Y_gen_mag, Y_gen_angle]=Ygen(Y_mat_ex,m,n)
for i=1:m
     for j=1:m
     Y11(i,j)=Y_mat_ex(i,j);
     end
 end
  for i=1:m
     for j=1:n
     Y12(i,j)=Y_mat_ex(i,j+4);
     end
  end
 for i=1:n
     for j=1:m
     Y21(i,j)=Y_mat_ex(i+4,j);
     end
 end
 for i=1:n
     for j=1:n
     Y22(i,j)=Y_mat_ex(i+4,j+4);
     end
 end
 
Y_gen=Y11-Y12*inv(Y22)*Y21; 
Y_gen_mag=abs(Y_gen);
Y_gen_angle=angle(Y_gen); %returns the phase angles in radians