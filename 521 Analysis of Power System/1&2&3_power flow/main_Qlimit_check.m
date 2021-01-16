clc; clear all
%% read grid data from file
file_name='ieee14cdf.txt';
[S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data,V_mag_final,V_ang_final,P_load,Q_load,P_gen,Q_gen,Qmax,Qmin]=read_data(file_name)    
%% power flow solution
% [PQ,PV,nPQ,nPV,P_cal, Q_cal,V_result]= LU_NRpowerflow(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data,V_mag_final,V_ang_final,P_load,Q_load,P_gen,Q_gen,Qmax,Qmin)
[PQ,PV,nPQ,nPV,P_cal, Q_cal,V_result,iter]= fast_decoupled_PF(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data,V_mag_final,V_ang_final,P_load,Q_load,P_gen,Q_gen,Qmax,Qmin)
 
%% check Qlimit for PV buses
Qmax(2)=25/S_Base; % set Qmax limit for bus 2 to 25Mvar, so that the Q_cal is greater than the Q limit
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
% [PQ,PV,nPQ,nPV,P_cal, Q_cal,V_result]= LU_NRpowerflow(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data,V_mag_final,V_ang_final,P_load,Q_load,P_gen,Q_gen,Qmax,Qmin)
[PQ,PV,nPQ,nPV,P_cal, Q_cal,V_result,iter]= fast_decoupled_PF(S_Base,No_of_Buses,No_of_Lines,Bus_data,Line_data,V_mag_final,V_ang_final,P_load,Q_load,P_gen,Q_gen,Qmax,Qmin)
 end

        



    


 
 
  
 
 
 
 
 
    
 
 
 
 
 
         
         
         






    
