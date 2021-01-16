function [Scheme2_order]=Scheme2(A)
% example for testing
% clc; clear all;
% 10 nodes example
% A=[1,1,0,1,0,0,0,1,0,0;1,1,0,0,0,0,1,0,0,1;0,0,1,1,1,0,1,1,0,1;1,0,1,1,1,0,0,0,0,0;0,0,1,1,1,1,1,0,0,1;0,0,0,0,1,1,1,0,0,0;0,1,1,0,1,1,1,1,1,0;1,0,1,0,0,0,1,1,0,0;0,0,0,0,0,0,1,0,1,0;0,1,1,0,1,0,0,0,0,1];
%%
n=length(A);
Scheme2_order=[];
index_fills=(1:n);
  %% calculate fills for each node if eliminated
while ~isempty(index_fills)
    m=length(index_fills);  %each elimination alternative is considered
    fills=zeros(m,1);
    for i=1:m
        ii=index_fills(i);
        A1_nnz=find(A(ii,:)~=0);
        if length(A1_nnz)>2 
            for p=1:length(A1_nnz)-1   
                for q=p+1:length(A1_nnz)
                    if (A1_nnz(p)~=ii)&&(A1_nnz(q)~=ii)&&A(A1_nnz(p), A1_nnz(q))==0  %build connection for other nodes
                    fills(i)=fills(i)+1;
                    end
                end
            end
        end
    end
    %% find the least fills
    k_fills=unique(fills); 
    index_least_fills=find(fills==k_fills(1));
    if length(index_least_fills)==1
        kk=index_least_fills(1);
        node_eliminated=index_fills(kk);
    else
       degrees=[];
       for i=1:length(index_least_fills)
           ii=index_least_fills(i);
           kk=index_fills(ii); % find node number
           degrees=[degrees,sum(A(:,kk)~=0)-1];
       end
       k=unique(degrees);
       index_least_degrees=find(degrees==k(1));       
       ii=index_least_degrees(1);
       jj=index_least_fills(ii);
       node_eliminated=index_fills(jj);
    end      
    Scheme2_order=[Scheme2_order,node_eliminated];
   
    %% edges created
    A2_nnz=find(A(node_eliminated,:)~=0);
    if length(A2_nnz)>2 
        for i=1:length(A2_nnz)-1  
            ii=node_eliminated;
            for j=i+1:length(A2_nnz)
                if (A2_nnz(i)~=ii)&&(A2_nnz(j)~=ii)&& A(A2_nnz(i), A2_nnz(j))==0  %build connection for other nodes
                    A(A2_nnz(i), A2_nnz(j))=1;
                    A(A2_nnz(j), A2_nnz(i))=1;
                end
            end
        end
    end
    A(node_eliminated,:)=0; %delete the node with the least degree
    A(:,node_eliminated)=0;
    k=find(index_fills==node_eliminated);
    index_fills(k)=[];       
end
end


        
            
            
            
            
            
        
        
    
    

