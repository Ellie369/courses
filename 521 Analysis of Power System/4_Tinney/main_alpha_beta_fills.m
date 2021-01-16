
% only for testing this script
clc; clear all;
% 10 nodes example
A=[11 12 0 14 0 0 0 18 0  0; %1
   21 22 0 0 0 0 27 0 0 210; %2
   0 0 33 34 35 0 37 38 0 310; %3 
   41 0 43 44 45 0 0 0 0  0; %4
   0 0 53 54 55 56 57 0 0 510; %5
   0 0 0 0 65 66 67 0 0  0; %6
   0 1 2 0 75 76 77 78 79  0; %7
   1 0 83 0 0 0 87 88 0  0; %8
   0 0 0 0 0 0 97 0 99  0; %9
   0 102 103 0 105 0 0 0 0 1010; %10
   ];
n=length(A);
%% before ordering
[Q_0, alpha_0,beta_0]=Qmatrix_alpha_beta(A)

% figure;
subplot(1,2,1);
spy(A);
title('A without ordering')

subplot(1,2,2);
spy(Q_0);
title('Q without ordering')
%% scheme0
[Scheme0_order]=Scheme0(A);
A_scheme0=[];
for i=1:n
    ii=Scheme0_order(i);
    for j=1:n
        jj=Scheme0_order(j)
        A_scheme0(i,j)=A(ii,jj);
    end
end
[Q_scheme0, alpha_scheme0,beta_scheme0]=Qmatrix_alpha_beta(A_scheme0)

% figure;
subplot(1,2,1);
spy(A_scheme0);
title('A with Scheme0')

subplot(1,2,2);
spy(Q_scheme0);
title('Q with Scheme0')
%% scheme1
[Scheme1_order]=Scheme1(A);
A_scheme1=[];
for i=1:n
    ii=Scheme1_order(i);
    for j=1:n
        jj=Scheme1_order(j)
        A_scheme1(i,j)=A(ii,jj);
    end
end
[Q_scheme1, alpha_scheme1,beta_scheme1]=Qmatrix_alpha_beta(A_scheme1)

% figure;
subplot(1,2,1);
spy(A_scheme1);
title('A with Scheme1')

subplot(1,2,2);
spy(Q_scheme1);
title('Q with Scheme1')

%% scheme2
[Scheme2_order]=Scheme2(A);
A_scheme2=[];
for i=1:n
    ii=Scheme2_order(i);
    for j=1:n
        jj=Scheme2_order(j)
        A_scheme2(i,j)=A(ii,jj);
    end
end
[Q_scheme2, alpha_scheme2,beta_scheme2]=Qmatrix_alpha_beta(A_scheme2)

figure;
subplot(1,2,1);
spy(A_scheme2);
title('A with Scheme2')

subplot(1,2,2);
spy(Q_scheme2);
title('Q with Scheme2')


