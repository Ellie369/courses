
% clc; clear all
% close all
function [PQ,nPQ,PV,nPV,Theta,Y_mag,V_mag,V_Delta,P_gen_cal,Q_gen_cal, V_result,P_load, Q_load]= NR_power_flow(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data)
% %% read grid data from file
% file_name='b_kundur_system.txt';
% [S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data]=read_data(file_name);

 V_mag_final=Bus_data(:,5); %magnitude of final voltage in p.u.
 V_ang_final=Bus_data(:,6); % change angle of final voltage from degrees to radians.
 P_load=Bus_data(:,7)/S_Base; %active power of Load in p.u.
 Q_load=Bus_data(:,8)/S_Base; %reactive power of Load in p.u.
 P_gen=Bus_data(:,9)/S_Base; %active power of generation in p.u.
 Q_gen=Bus_data(:,10)/S_Base; %reactive power of genertion in p.u.

 Qmax = Bus_data(:,13)/S_Base;     % Maximum Reactive Power Limit
 Qmin = Bus_data(:,14)/S_Base;      % Minimum Reactive Power Limit

 %% type of buses
 slack=find(Bus_data(:,4)==3); %slack bus
 PQ=find(Bus_data(:,4)==0|Bus_data (:,4) == 1); %PQ bus
 PV=find(Bus_data(:,4)==2); %PV bus
 nslack=length(slack);
 nPQ=length(PQ);
 nPV=length(PV);

 
  %% form Y_matrix
 [Y_mat,Theta,Y_mag,B,G]=y_bus(Bus_data,Line_data,No_of_Buses,No_of_Lines); %form Y_matrix
 
 % load data
 
 %% initialize parameters
% V_mag=Bus_data(:,5);
% V_Delta=Bus_data(:,6)./180*pi; %initial angle

V_mag=Bus_data(:,12); %initial voltage
V_mag(~V_mag) = 1; % replace all 0 with 1 for voltage magnitude
V_Delta=zeros(No_of_Buses,1); %initial angle

% % initialize ZIP lodad model
% P_p=0.5; P_q=0.5; %constant power percent
% I_p=0; I_q=0;   %constant current percent
% Z_p=0.5; Z_q=0.5; % constant impedance percent
% 
% PL0=P_p* P_load;
% M0=(I_p* P_load)./V_mag;
% G0=(Z_p*P_load)./(V_mag.^2);
% QL0=P_q* Q_load;
% H0=(I_q* Q_load)./V_mag; 
% B0=(Z_q*Q_load)./(V_mag.^2);

P_sch=P_gen-P_load; %net power scheduled at a bus
Q_sch=Q_gen-Q_load;
 dif_Voltage=zeros(No_of_Buses-1+nPQ,1);
 
 %% Newton-Raphson calculation start
 
 tol_max=0.01; % iteration tolerance of the solution
 iter=0; %times of iteration
 tol=1; %initial value of tol
 
 while (tol>tol_max)
%      P_load=PL0+M0.*V_mag + G0.*(V_mag.^2);
%      Q_load=QL0+H0.*V_mag + B0.*(V_mag.^2);
%      P_sch=P_gen-P_load; %net power scheduled at a bus
%      Q_sch=Q_gen-Q_load;
     
    [P_cal,Q_cal]=cal_PQ(V_mag,Y_mag,Theta,V_Delta,No_of_Buses); %calculate the power at buses     
    [dif_PQ]=difference_PQ(P_sch,Q_sch,P_cal,Q_cal,PQ,nPQ); % mismatches vector
    [J]=Jacobian_matrix(V_mag,P_cal,Q_cal,Y_mag,Theta, V_Delta,No_of_Buses,PQ,nPQ,B,G); %call the Jacobian matrix
    dif_Voltage=inv(J)*dif_PQ; %get correction vector
    dif_D=dif_Voltage(1:No_of_Buses-1); % angle correction vector
    dif_V=dif_Voltage(No_of_Buses:end); % magnitude correction vector
    
    V_Delta(2:end)=V_Delta(2:end)+dif_D; %correct the results, angle and voltage
    V_mag(PQ)=V_mag(PQ)+dif_V;
    tol=max(abs(dif_PQ));
    iter=iter+1;
    
    
 
 end
    if iter>12
        disp('bad, not converge');
    else
       
        V_result=[V_mag,rad2deg(V_Delta)];
        P_gen_cal=P_load+P_cal;
        Q_gen_cal=Q_load+Q_cal;
            
    end  
% %% verify the calculation results
% if max(V_result-[V_mag_final,V_ang_final])<=0.1
%     fprintf('Congratulation,converge!, times of iteration=%d.\n',iter);
% else
%     disp('converge, but results are not correct, please go back to check!');
% end
    
end    


 
 
  
 
 
 
 
 
    
 
 
 
 
 
         
         
         






    
