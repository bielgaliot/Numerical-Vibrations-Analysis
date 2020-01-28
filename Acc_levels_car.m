clc, clear all, close all;
tic

m=480; % [kg]
J=410; % [kgm^2] 
m1=53; % [kg]
m2=53; % [kg]

l1=1; % [m]
l2=0.8; % [m]

k=100000; % [N/m] 
kr=100000; % [N/m] 
kt=80000; % [N/m]

c=5000; % [Ns/m] 
ct=4000; % [Ns/m] 
cr=5000; % [Ns/m] 
d=1.2; % [m] 


M=diag([m J m1 m2]);

K=[2*k -k*l1+k*l2 -k -k;...
   -k*l1+k*l2 k*l1^2+k*l2^2+kr k*l1 -k*l2;...
   -k k*l1 k+kt 0;...
   -k -k*l2 0 kt+k];

C=[2*c -c*l1+c*l2 -c -c;...
   -c*l1+c*l2 c*l1^2+c*l2^2+cr c*l1 -c*l2;...
   -c c*l1 c+ct 0;...
   -c -c*l2 0 ct+c];

N=length(M); % Number of DOFs of the problem.
%% Normal modes
DOF=length(M);
[Phi, Lambda]=eig(K,M)
fn=sqrt(Lambda)/(2*pi)
for i=1:DOF
    Phi_norm(:,i)=Phi(:,i)/Phi(1,i); 
end
Phi_norm
Phi

%% Acc Analisis modal
we=100;

F=[-800 -800*l2 0 0].';
F_mod=Phi.'*F;
K_mod= Phi.'*K*Phi
C_mod= Phi.'*C*Phi;
C_mod_diag=diag(diag(C_mod)) 
M_mod= Phi.'*M*Phi
      
    Z_mod=(-we^2*M_mod+1i*we*C_mod_diag+K_mod);
    Q=Z_mod\F_mod;
    
X_modal_meth=Phi*Q;
Xacc_modal=we^2*X_modal_meth;

Xacc_modal_A=abs(Xacc_modal(1)-d*Xacc_modal(2))
Xacc_modal_B=abs(Xacc_modal(1)+d*Xacc_modal(2))

%% Acc Analisis direct

Z=(-we^2*M+1i*we*C+K);
X=Z\F;

Xacc_directo=-we^2*X;

Xacc_directo_A=abs(Xacc_directo(1)-d*Xacc_directo(2))
Xacc_directo_B=abs(Xacc_directo(1)+d*Xacc_directo(2))

%% Acc Analisis direct with base excitation

wb=20;
F=[0 0 0.2*(kt+1i*wb*ct) 0.1*(kt+1i*wb*ct)].';
Z=(-wb^2*M+1i*wb*C+K);
X=Z\F;

Xacc_directo_Y=-wb^2*X

Xacc_directo_A_Y=abs(Xacc_directo_Y(1)-d*Xacc_directo_Y(2))
Xacc_directo_B_Y=abs(Xacc_directo_Y(1)+d*Xacc_directo_Y(2))

%% Acc in a and b with c=0 for range of w

w=linspace(0, 70, 1000);



for n=1:length(w)
  
    F=[0 0 0.2*(kt) 0.1*kt].';
    
    Z=-w(n)^2*M+K;
    X=Z\F;
    
    X_acc=-w(n)^2*X;
    X_acc_A(n)=abs(X_acc(1)-d*X_acc(2));
    X_acc_B(n)=abs(X_acc(1)+d*X_acc(2));
    
  
end
semilogy(w,X_acc_A)%DoF #1
hold on
grid
semilogy(w,X_acc_B)
hold on










