
clc; 
clear all
%% read grid data from file and power flow solution
file_name='ieee14cdf.txt';
[S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data]=read_data(file_name);    
[PQ,PV,nPQ,nPV,V_mag,V_Delta,P_cal,Q_cal,Y_mat,Theta,Y_mag,B,G,Pij,Qij,Pji,Qji]= LU_NRpowerflow(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data);
% True_value=[V_mag;P_cal;Q_cal;Pij;Pji;Qij;Qji]; 

%% add noise to measurements
z = [V_mag;P_cal;Q_cal;Pij;Pji;Qij;Qji]; 
V_noise = 0.01;
bus_noise = 0.015;
Line_noise = 0.02;
noise_Vm = normrnd(0,V_noise,length(V_mag),1);
noise_Inj = normrnd(0,bus_noise,length([P_cal;Q_cal]),1);
noise_flow = normrnd(0,Line_noise,length([Pij;Pji;Qij;Qji]),1);
z = z + [noise_Vm;noise_Inj;noise_flow];

z(84) = z(84) + 0.2;
% z(85) = z(85) + 0.2;
% z(86) = z(86) + 0.2;

W = diag([1/(V_noise^2)*ones(length(V_mag),1);1/(bus_noise^2)*ones(length([P_cal;Q_cal]),1);1/(Line_noise^2)*ones(length([Pij;Pji;Qij;Qji]),1)]); % weighting matrix
iter = 0;
tol = 1;

while(tol > 1e-3)
   [h]=h_updata(V_mag,V_Delta,Line_data,No_of_Buses,No_of_Lines,Theta,Y_mag);
    e = z - h;
   H=H_Matrix(V_mag,V_Delta,Line_data,No_of_Buses,No_of_Lines,Theta,Y_mag,B,G);
    Gm = H'*W*H;
    X_dif = inv(Gm)*(H'*W*e);
    V_Delta(2:end) = V_Delta(2:end) + X_dif(1:No_of_Buses-1);
    V_mag = V_mag + X_dif(No_of_Buses:end);
    iter = iter + 1;
    tol = max(abs(X_dif));
    fprintf('State Estimator: The maximum update of state variables in %d iteration is %f\n',iter,tol);
end
fprintf('State Estimator: Converged in %d iterations!\n',iter);
result = [V_Delta(2:end);V_mag];
h=h_updata(V_mag,V_Delta,Line_data,No_of_Buses,No_of_Lines,Theta,Y_mag);

%% bad data detection 
e = z - h;
f_hat = e'*W*e;
fprintf(' f for this state estimator is %d\n',f_hat);
if f_hat > 129.973
    H=H_Matrix(V_mag,V_Delta,Line_data,No_of_Buses,No_of_Lines,Theta,Y_mag,B,G);
    Gm = H'*W*H;
    R = diag(1./diag(W));
    R_prime = diag(R - H*(Gm\H'));
    Var = e./sqrt(R_prime);
    plot(Var,'.-');
    Var=abs(Var);
    bad_data=find(Var==max(Var));
    fprintf('Detected bad data is from z(%d)\n',bad_data);
end
% %% bad data detection for sevetal measurements
% e = z - h;
% f_hat = e'*W*e;
% fprintf(' f for this state estimator is %d\n',f_hat);
% while f_hat > 129.973
%     H=H_Matrix(V_mag,V_Delta,Line_data,No_of_Buses,No_of_Lines,Theta,Y_mag,B,G);
%     Gm = H'*W*H;
%     R = diag(1./diag(W));
%     R_prime = diag(R - H*(Gm\H'));
%     Var = e./sqrt(R_prime);
%     plot(Var,'.-');
%     Var=abs(Var);
%     bad_data=find(Var==max(Var));
%     fprintf('Detected bad data is from z(%d)\n',bad_data);
%     
% end
% 

    

%     e_2=e; e_2(bad_data)=[];
%    W_2=W; W_2(bad_data,:)=[];W_2(:,bad_data)=[];
%     f_hat = e_2'*W_2*e_2;
% else
%     break
% end