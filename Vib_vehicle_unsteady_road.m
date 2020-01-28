clc
clear all
close all


% Mechanical or geometrical parameters of the model
m1=75; % [kg]
m2=75; % [kg]
m3=900; % [kg]
J=1900; % [kgm?2]
k=2e4; % [N/m]
k_wp=8e5; % [N/m] 
c=3500; % [Ns/m] 
l=1; % [m]
l3=0.3; % [m]
%% QUESTION 1. Mass, damping and stiffness matrices. 
% Mass, damping and stiffness matrices.
M=diag([m3 J m1 m2]);

K=[2*k,0,-k,-k;...
   0 ,2*k*l^2, k*l,-k*l;...
   -k, k*l ,k_wp+k, 0;...
   -k,-k*l,0,k_wp+k];

C=[2*c,0,-c,-c;...
    0,2*c*l^2,c*l,-c*l;...
    -c,c*l,c,0;...
    -c,-c*l,0,c]; 

N=length(M); % DOF
%% Normal modes

%Contact
DOF=length(M);
[Phi, Lambda]=eig(K,M);
fn=sqrt(Lambda)/(2*pi)
for i=1:DOF
    Phi_norm(:,i)=Phi(:,i)/Phi(1,i); 
end
Phi_norm
Phi*10^3






%Not in contact
Knc=[2*k,0,-k,-k;...
   0 ,2*k*l^2, k*l,-k*l;...
   -k, k*l,k, 0;...
   -k,-k*l,0,k];

DOF=length(M);
[PhiNC, LambdaNC]=eig(Knc,M);
fnNC=sqrt(LambdaNC)/(2*pi);
for i=1:DOF
    Phi_norm(:,i)=Phi(:,i)/PhiNC(1,i); 
end
Phi_norm;
PhiNC*10^3


%% QUESTION 3. Normal modes.
F=1000; % [N]

we3=2500*(2*pi)/60; % [rad/s]

F3=[F F*l3 0 0].'; % [N]

% Direct Method

X3=(K+1i*we3*C-we3^2*M)\F3;

% Vibration velocity related to the motion of the vehicle over wheel 1.
X3_vehc=X3(1)-X3(2)*l;
abs(X3_vehc);

X3_veldB_vehc=20*log10(abs(-1i*we3*X3_vehc)/1e-9)

%% QUESTION 3. Excitació paviment

%Mètode directe;
 Xp=0.1;
 lambda=5;
 F=[0 0 k_wp*Xp k_wp*Xp*exp(1i*2*pi*(2*l)/lambda)].';
    
 v=90/3.6;
 we=10*pi;
 X=(K+1i*we*C-we^2*M)\F;

 X_despl=abs(X(1)) %I[m]

%% Are dampers working properly?

w=linspace(0, 2*pi*40, 1000);
F=[0 0 k_wp*Xp k_wp*Xp*exp(1i*2*pi*(2*l)/lambda)].'


X_dir=zeros(length(w),1);

for n=1:length(w)
  for j=1:4
      
    Z=(-w(n)^2*M+1i*w(n)*C+K);
    X=Z\F;
    X_dir(j,n)=abs(X(j));

  end
end

plot(w/(2*pi),X_dir(1,:))%modeshape 1, 
hold on

grid
 
 
 




