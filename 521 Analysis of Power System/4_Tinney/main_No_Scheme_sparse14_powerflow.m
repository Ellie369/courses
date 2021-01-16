% function [nPV, PV, P_cal, Q_cal,V_result]= LU_NRpowerflow(Bus_data, Line_data,S_Base,No_of_Buses,No_of_Lines,P_gen,P_load,Q_gen,Q_load,V_mag_final,V_ang_final)
%% Jacobian matrix is stored using sparse technique 
%% read data
clc; clear all
file_name='ieee14cdf.txt';
[S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data,]=read_data(file_name); 
%% analyze data  
 V_mag_final=Bus_data(:,5); %magnitude of final voltage in p.u.
 V_ang_final=Bus_data(:,6); % change angle of final voltage from degrees to radians.
 P_load=Bus_data(:,7)/S_Base; %active power of Load in p.u.
 Q_load=Bus_data(:,8)/S_Base; %reactive power of Load in p.u.
 P_gen=Bus_data(:,9)/S_Base; %active power of generation in p.u.
 Q_gen=Bus_data(:,10)/S_Base; %reactive power of genertion in p.u.
 Qmax = Bus_data(:,13)/S_Base;     % Maximum Reactive Power Limit
 Qmin = Bus_data(:,14)/S_Base;     % Minimum Reactive Power Limit
 P_sch=P_gen-P_load; %net power scheduled at a bus
 Q_sch=Q_gen-Q_load;
%% type of buses
 slack=find(Bus_data(:,4)==3); %slack bus
 PQ=find(Bus_data(:,4)==0|Bus_data (:,4) == 1); %PQ bus
 PV=find(Bus_data(:,4)==2); %PV bus
 nslack=length(slack);
 nPQ=length(PQ);
 nPV=length(PV);

  %% form Y_matrix
 [Y_mat,Theta,Y_mag,B,G]=y_bus(Bus_data,Line_data,No_of_Buses,No_of_Lines); %form Y_matrix
 
 %% initialize parameters
 
V_mag=Bus_data(:,12); %initial voltage
V_mag(~V_mag) = 1; % replace all 0 with 1 for voltage magnitude
V_Delta=zeros(No_of_Buses,1); %initial angle
dif_Voltage=zeros(No_of_Buses-1+nPQ,1);
 
 %% Newton-Raphson calculation start
 
 tol_max=0.01; % iteration tolerance of the solution
 iter=0; %times of iteration
 tol=1; %initial value of tol
 
 while (tol>tol_max)
    [P_cal,Q_cal]=cal_PQ(V_mag,Y_mag,Theta,V_Delta,No_of_Buses); %calculate the power at buses     
    [dif_PQ]=difference_PQ(P_sch,Q_sch,P_cal,Q_cal,PQ,nPQ); % mismatches vector
    [J]=Jacobian_matrix(V_mag,P_cal,Q_cal,Y_mag,Theta, V_Delta,No_of_Buses,PQ,nPQ,B,G); %call the Jacobian matrix   
    %dif_Voltage=inv(J)*dif_PQ; %get correction vector
%     spy(J);
    %% Tinney 0/1/2, so as to get the ordering of nodes
%     [Scheme0_order]=Scheme0(J);
%     [Scheme1_order]=Scheme1(J);
%     [Scheme2_order]=Scheme2(J);

%% Sparse storage of Jacobian matrix
%     [J_ordered,FIC,FIR,NCOL,NROW,NIC,NIR,dif_PQ_ordered]=Sparse_storage(Scheme0_order,J,dif_PQ);
%     [J_ordered,FIC,FIR,NCOL,NROW,NIC,NIR,dif_PQ_ordered]=Sparse_storage(Scheme1_order,J,dif_PQ);
%     [J_ordered,FIC,FIR,NCOL,NROW,NIC,NIR,dif_PQ_ordered]=Sparse_storage(Scheme2_order,J,dif_PQ);

    %% LU factorization 

%     [dif_Voltage_ordered]=LU_sparse(J_ordered,FIC,FIR,NCOL,NROW,NIC,NIR,dif_PQ_ordered); 
[dif_Voltage]=LU_factor_PQ(J,dif_PQ);
%% get voltage vector in the original sequence

%     dif_Voltage(Scheme0_order)=dif_Voltage_ordered;
%     dif_Voltage(Scheme1_order)=dif_Voltage_ordered;
%     dif_Voltage(Scheme2_order)=dif_Voltage_ordered;
    
    dif_D=dif_Voltage(1:No_of_Buses-1); % angle correction vector
    dif_V=dif_Voltage(No_of_Buses:end); % magnitude correction vector
    
    V_Delta(2:end)=V_Delta(2:end)+dif_D; %correct the results, angle and voltage
    V_mag(PQ)=V_mag(PQ)+dif_V;
    tol=max(abs(dif_PQ));
    iter=iter+1;
 end
    if iter>5
        disp('bad, not converge');
    else
        V_result=[V_mag,rad2deg(V_Delta)];        
    end  
%% verify the calculation results
if max(V_result-[V_mag_final,V_ang_final])<=0.1
    fprintf('Congratulation,converge!, times of iteration=%d.\n',iter);
else
    disp('converge, but results are not correct, please go back to check!');
end
    
%% calculate alpha beta and fills without Scheme
[Q_J, alpha,beta]=Qmatrix_alpha_beta(J);
figure;
subplot(1,2,1);
spy( J);
title('J without Scheme')

subplot(1,2,2);
spy(Q_J);
title('Q without Scheme')
% %%  Scheme0/1/2, calculate alpha beta and fills
% n=length(J);
% J_ordered_matrix=[];
% for i=1:n
%     ii=Scheme0_order(i);
%     for j=1:n
%         jj=Scheme0_order(j);
%         J_ordered_matrix(i,j)=J(ii,jj);
%     end
% end
% 
% [Q_J_ordered, alpha,beta]=Qmatrix_alpha_beta(J_ordered_matrix);
% alpha
% beta
% figure;
% subplot(1,2,1);
% spy( J_ordered_matrix);
% title('J with Scheme 0')
% 
% subplot(1,2,2);
% spy(Q_J_ordered);
% title('Q with Scheme 0')
 
 
  
 
 
 
 
 
    
 
 
 
 
 
         
         
         






    
