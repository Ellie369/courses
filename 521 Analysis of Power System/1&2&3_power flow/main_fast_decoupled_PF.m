function [PQ,PV,nPQ,nPV,P_cal, Q_cal,V_result,iter]= fast_decoupled_PF(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data,V_mag_final,V_ang_final,P_load,Q_load,P_gen,Q_gen,Qmax,Qmin,slack,PQ,PV,nPQ,nPV)
% read grid data from file---for independent run
% clc; clear all;
% file_name='ieee14cdf.txt';
% [S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data]=read_data(file_name);
%  V_mag_final=Bus_data(:,5); %magnitude of final voltage in p.u.
%  V_ang_final=Bus_data(:,6); % change angle of final voltage from degrees to radians.
%  P_load=Bus_data(:,7)/S_Base; %active power of Load in p.u.
%  Q_load=Bus_data(:,8)/S_Base; %reactive power of Load in p.u.
%  P_gen=Bus_data(:,9)/S_Base; %active power of generation in p.u.
%  Q_gen=Bus_data(:,10)/S_Base; %reactive power of genertion in p.u.
%  Qmax = Bus_data(:,13)/S_Base;     % Maximum Reactive Power Limit
%  Qmin = Bus_data(:,14)/S_Base;     % Minimum Reactive Power Limit
 %% type of buses
 slack=find(Bus_data(:,4)==3); %slack bus
 PV=find(Bus_data(:,4)==2); %PV bus index
 PQ=find(Bus_data(:,4)==0|Bus_data (:,4) == 1); %PQ bus index
 nslack=length(slack);
 nPQ=length(PQ);
 nPV=length(PV);

  %% form Y_matrix
 [Y_mat,Theta,Y_mag,B,G]=y_bus(Bus_data,Line_data,No_of_Buses,No_of_Lines); %form Y_matrix
 
%% from Y_matrix obtain B_mat_dec that used in decouple power flow
[B1,B2]=B_mat(B,nPQ,PQ,No_of_Buses);

%% initialize parameters
V_mag=Bus_data(:,12); %initial voltage
V_mag(~V_mag) = 1; % replace all 0 with 1 for voltage magnitude
V_Delta=zeros(No_of_Buses,1); %initial angle
 P_sch=P_gen-P_load; %net power scheduled at a bus
 Q_sch=Q_gen-Q_load;
 dif_Voltage=zeros(No_of_Buses-1+nPQ,1);
 
 %% Decoupled power flow calculation start
 
 tol_max=0.01; % iteration tolerance of the solution
 iter=0; %times of iteration
 tol=1; %initial value of tol
 
 while (tol>tol_max)
    [P_cal,Q_cal]=cal_PQ(V_mag,Y_mag,Theta,V_Delta,No_of_Buses); %calculate the power at buses       
    % calculate ?P/V
    dif_P_all=P_sch-P_cal;
    dif_P=dif_P_all(2:end);
    for i=2:No_of_Buses
        dif_P_V(i-1)=dif_P_all(i)/V_mag(i); 
    end
    % calculate voltage angle correction vector using LU factorization
    [dif_D]=LU_factor_PQ(B1,dif_P_V); 
    
    % calculate ?Q/V
    dif_Q_all=Q_sch-Q_cal;
    dif_Q=zeros(nPQ,1);
    dif_Q_V=zeros(nPQ,1);
    for i=1:nPQ
        dif_Q(i)=dif_Q_all(PQ(i));
        dif_Q_V(i)=dif_Q_all(PQ(i))/V_mag(PQ(i));
    end
    % calculate voltage magnitude correction vector using LU factorization
     [dif_V]=LU_factor_PQ(B2,dif_Q_V); 
     
     %correct the results, voltage magnitude and angle
     V_mag(PQ)=V_mag(PQ)+dif_V;
     V_Delta(2:end)=V_Delta(2:end)+dif_D;
     dif_PQ=[dif_P;dif_Q];
     tol=max(abs(dif_PQ));
     iter=iter+1;     
 end
     if iter>10
        disp('bad, not converge');
    else
        V_result=[V_mag,rad2deg(V_Delta)];        
    end  
%% verify the calculation results
if max(V_result-[V_mag_final,V_ang_final])<=0.5
    fprintf('Congratulation,converge!, times of iteration=%d.\n',iter);
else
    disp('converge, but results are not correct, please go back to check!');
end
end    
    


 
 
  
 
 
 
 
 
    
 
 
 
 
 
         
         
         






    
