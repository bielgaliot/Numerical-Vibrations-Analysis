%% Problem 1


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





%% Due to rotative force with damping

wr=linspace(0, 2*pi*50, 1000);

K_mod= Phi.'*K*Phi;
C_mod= Phi.'*C*Phi;
C_mod_diag=diag(diag(C_mod)); 
M_mod= Phi.'*M*Phi;

F=[-1i -1i*l2 0 0].';
F_mod=Phi.'*F;

Q_mod=zeros(length(wr),1);


for n=1:length(wr)
  for j=1:4 
    Z_mod=(-wr(n)^2*M_mod+1i*wr(n)*C_mod_diag+K_mod);
    Q=Z_mod\F_mod;
    Q_mod(j,n)=abs(Q(j)); 

  end
end

figure
semilogy(wr/(2*pi),Q_mod(1,:)) %DoF #1
hold on
semilogy(wr/(2*pi),Q_mod(2,:))
hold on
semilogy(wr/(2*pi),Q_mod(3,:))
hold on
semilogy(wr/(2*pi),Q_mod(4,:))
legend('MS1','MS2','MS3', 'MS4')
xlabel 'Frequency (Hz)'
ylabel '|Q|'
title 'Rotative Force Damped'
grid

%EXPLANATIONS WRITTEN BY HAND IN THE PIECES OF PAPER HANDED IN WITH THE
%EXAM
%% Due to rotative force without damping

wr=linspace(0, 2*pi*50, 1000);

K_mod= Phi.'*K*Phi;
C_mod= Phi.'*C*Phi;
C_mod_diag=diag(diag(C_mod)*0); %Damping is made 0 by multiplying by 0.
M_mod= Phi.'*M*Phi;

F=[-1i -1i*l2 0 0].';
F_mod=Phi.'*F;
    
Q_mod=zeros(length(wr),1);


for n=1:length(wr)
  for j=1:4 
    Z_mod=(-wr(n)^2*M_mod+1i*wr(n)*C_mod_diag+K_mod);
    Q=Z_mod\F_mod;
    Q_mod(j,n)=abs(Q(j)); 

  end
end

Q_mod(j,n)


figure
semilogy(wr/(2*pi),Q_mod(1,:)) %DoF #1
hold on
semilogy(wr/(2*pi),Q_mod(2,:))
hold on
semilogy(wr/(2*pi),Q_mod(3,:))
hold on
semilogy(wr/(2*pi),Q_mod(4,:))
legend('MS1','MS2','MS3', 'MS4')
xlabel 'Frequency (Hz)'
ylabel '|Q|'
grid
title 'Rotative Force Undamped'

%% Due to base excitation with damping

%WITH DAMPING

wb=linspace(0, 2*pi*50, 1000);

K_mod= Phi.'*K*Phi;
C_mod= Phi.'*C*Phi;
C_mod_diag=diag(diag(C_mod)); 
M_mod= Phi.'*M*Phi;

Q_mod=zeros(length(wb),1);


for n=1:length(wb)
  for j=1:4 
      
    F=[0 0 0.2*(kt+1i*wb(n)*ct) 0.1*(kt+1i*wb(n)*ct)].';
    F_mod=Phi.'*F;
    Z_mod=(-wb(n)^2*M_mod+1i*wb(n)*C_mod_diag+K_mod);
    Q=Z_mod\F_mod;
    Q_mod(j,n)=abs(Q(j)); 

  end
end

%Q_mod(j,n)

figure
semilogy(wb/(2*pi),Q_mod(1,:)) %DoF #1
hold on
semilogy(wb/(2*pi),Q_mod(2,:))
hold on
semilogy(wb/(2*pi),Q_mod(3,:))
hold on
semilogy(wb/(2*pi),Q_mod(4,:))
legend('MS1','MS2','MS3', 'MS4')
xlabel 'Frequency (Hz)'
ylabel '|Q|'
title 'Base Excitation Damped'
grid

%% Due to base excitation without damping

%WITHOUT DAMPING

wb=linspace(0, 2*pi*50, 1000);

K_mod= Phi.'*K*Phi;
C_mod= Phi.'*C*Phi;
C_mod_diag=diag(diag(C_mod)*0); %Damping is made 0 by multiplying by 0.
M_mod= Phi.'*M*Phi;

Q_mod=zeros(length(wb),1);


for n=1:length(wb)
  for j=1:4 
      
    F=[0 0 0.2*(kt+1i*wb(n)*ct) 0.1*(kt+1i*wb(n)*ct)].';
    F_mod=Phi.'*F;
    Z_mod=(-wb(n)^2*M_mod+1i*wb(n)*C_mod_diag+K_mod);
    Q=Z_mod\F_mod;
    Q_mod(j,n)=abs(Q(j)); 

  end
end

%Q_mod(j,n)

figure
semilogy(wb/(2*pi),Q_mod(1,:)) %DoF #1
hold on
semilogy(wb/(2*pi),Q_mod(2,:))
hold on
semilogy(wb/(2*pi),Q_mod(3,:))
hold on
semilogy(wb/(2*pi),Q_mod(4,:))
legend('MS1','MS2','MS3', 'MS4')
xlabel 'Frequency (Hz)'
ylabel '|Q|'
title 'Base Excitation Undamped'
grid


%% Problem 2
m=1000; % [kg]
J=810; % [kgm^2] 

l1=1; % [m]
l2=1.5; % [m]
l3=0.5; % [m]

kf=18e3; % [N/m]
kr=22e3; % [N/m] 


M=diag([m J]); 

K=[kf+kr kr*l2-kf*l1;...
   kr*l2-kf*l1 kf*l1^2+kr*l2^2];

C=[0 0;0 0];


N=length(M); % Number of DOFs of the problem.

%% Normal modes
%A=M\K;

[Phi2, Lambda2]=eig(K,M)
fn=sqrt(Lambda2)/(2*pi)

%% Acc en punto a en función de la velocidad
lambda=10;
A=0.05; % [m]

F=[kf*A+kr*A*exp(1i*2*pi*(l1+l2)/lambda) -kf*l1*A+kr*l2*A*exp(1i*2*pi*(l1+l2)/lambda)].';
    

v=linspace(20/3.6, 80/3.6, 1000);


for i=1:length(v)
   
        we=2*pi*(v(i)/lambda); 
        X=(-we^2*M+K)\F;
        X_A=abs((X(1)+(l2+l3)*X(2)));
        X_acc_A(i)=-(((v(i)/lambda)*2*pi)^2)*X_A;
   

end

figure
semilogy(v*3.6,abs(X_acc_A))
xlabel 'Car speed (km/h)'
ylabel('Acc at point A (m/s^2)')
title 'Acc at point A as a function of the velocity'
grid

%Big spikes at 54 and 33 km/h, the velocities which produce a frequency of
%excitation that corresponds with the natural frequencies of the system, calculated before in line 245. 
