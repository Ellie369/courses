%% make fast_decouple power flow calculation a function so as to check Q limit after power flow calculation
clc; clear all
%% read grid data from file
file_name='ieee14cdf.txt';
[S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data]=read_data(file_name);
 
[nPV, PV, P_cal, Q_cal,V_result]= LU_NRpowerflow(Bus_data, Line_data,S_Base,No_of_Buses,No_of_Lines,P_gen,P_load,Q_gen,Q_load,V_mag_final,V_ang_final);
 %% check Qlimit for PV buses
Q_cal(2)=80/S_Base; % set Q_cal for bus 2 to 80Mvar, so that the Q_cal is greater than the Q limit
while max(Q_cal(PV)>Qmax(PV))|| max(Q_cal(PV)<Qmin(PV))
    for i=1:nPV
        if Q_cal(PV(i))>Qmax(PV(i))
            Q_cal(PV(i))=Qmax(PV(i));
            Bus_data(PV(i),4)=0;
        elseif Q_cal(PV(i))<Qmin(PV(i))
            Q_cal(PV(i))=Qmin(PV(i));
            Bus_data(PV(i),4)=0;
        end
    end
[PQ,nPQ,nPV, PV, P_cal, Q_cal,V_result,iter]= fast_decoupled_PF(Bus_data, Line_data,S_Base,No_of_Buses,No_of_Lines,P_gen,P_load,Q_gen,Q_load,V_mag_final,V_ang_final);
 end

        



    


 
 
  
 
 
 
 
 
    
 
 
 
 
 
         
         
         






    
