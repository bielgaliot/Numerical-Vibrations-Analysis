clc, clear all, close all;
tic
mp1=4.8; % [kg]
mp2=7.2; % [kg]
mp3=6.4; % [kg]
Jp=1000; % [kg^2] 
kp1=18520; % [N/m] 
kp2=13912; % [N/m] 
kp3=4700; % [N/m] 
kpc=1500; % [N/m] 
cp1=970; % [Ns/m] 
cp2=750; % [Ns/m] 
cp3=250; % [Ns/m] 
a=0.5; % [m] 
zH=0.2; % [m]

M=diag([mp1 mp2 Jp mp3]);

K=[kp1+kp2 -kp2 0 0;...
    -kp2 kp2+kp3 -kp3*zH -kp3;...
    0 -kp3*zH kp2*a^2+kp3*zH^2 kp3*zH;...
    0 -kp3 kp3*zH kp3+kpc];

C=[cp1+cp2 -cp2 0 0;...
    -cp2 cp2+cp3 -cp3*zH -cp3;...
    0 -cp3*zH cp2*a^2+cp3*zH^2 cp3*zH;...
    0 -cp3 cp3*zH cp3];

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
%% Excitation
w=linspace(0, 2*pi*10, 1000);
F=[0 0 0 1].';
M_mod=Phi.'*M;
F_mod=Phi.'*F;
K_mod= Phi.'*K*Phi;
C_mod= Phi.'*C*Phi;
C_mod_diag=diag(diag(C_mod)); %multiply by 0 to remove damping without errors

Q_mod=zeros(length(w),1);

for n=1:length(w)
  for j=1:4 %Adaptar a DoF
      
    Z_mod=(-w(n)^2*M_mod+1i*w(n)*C_mod_diag+K_mod);
    Q=Z_mod\F_mod;
    Q_mod(j,n)=abs(Q(j));

  end
end
semilogy(w/(2*pi),Q_mod(1,:))%DoF #1
hold on
semilogy(w/(2*pi),Q_mod(2,:))
hold on
semilogy(w/(2*pi),Q_mod(3,:))
hold on
semilogy(w/(2*pi),Q_mod(4,:))
legend('MS1','MS2','MS3', 'MS4')
xlim([0 10])
grid
%ylim([])

%% Catenaria 2:

%Mètode directe;

 F_cat=[0 0 0 kpc*0.1].';
 v=linspace(0, 30, 1000);
 w=linspace(0, 30, 1000);
 X_save=zeros(4, length(v));
 lambda=16;

toc
Fpc=zeros(length(v),1); 

for i=1:length(v)
    
we=2*pi*(v(i)/lambda); 
X=(K+1i*we*C-we^2*M)\F_cat;
Fpc(i)=kpc*(X(4)-0.1);

end
figure
plot(v,abs(Fpc))
xlabel 'Velocitat (m/s)'
ylabel('F contacte (N)')
title 'Contact force as a function of speed'
grid

toc