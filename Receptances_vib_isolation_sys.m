clc, clear all, close all;
tic

m=20; % [kg]

k=4e6; % [N/m] 
 
c=3e3; % [Ns/m] 


M=diag([m 2*m]);

K=[2*k -k;...
   -k 2*k];
 
C=[2*c -c;...
   -c 2*c];

N=length(M); % Number of DOFsFor 
of the problem.

%% Normal modes
DOF=length(M);
[Phi, Lambda]=eig(K,M);
fn=sqrt(Lambda)/(2*pi)
for i=1:DOF
    Phi_norm(:,i)=Phi(:,i)/Phi(2,i); 
end
Phi_norm
Phi

%% Matriu de receptancies metode directe

%Receptance matrices with their plots as a function of freq
%Direct method
w=linspace(0, 2*pi*200, 1000);
H=zeros(N,N,length(w));

for n=1:length(w)
    Z=(-w(n)^2*M+1i*w(n)*C+K);
    H(:,:,n)=inv(Z);    
end

for j=1:length(w)
    
H11(j)=abs(H(1,1,j));
H12(j)=abs(H(1,2,j));
H21(j)=abs(H(2,1,j));
H22(j)=abs(H(2,2,j));

end

subplot(2,2,1)
plot(w/(2*pi), H11) 
ylabel('|H_{11}(w)|')
xlabel('Frequency (Hz)')
title 'Direct Receptance H_{11}'
grid

subplot(2,2,2)
plot(w/(2*pi), H12) %response in 1 due to unit force at 2
xlabel('Frequency (Hz)')
ylabel('|H_{12}(w)|')
title 'Direct Receptance H_{12}'
grid

subplot(2,2,3)
plot(w/(2*pi), H21)
xlabel('Frequency (Hz)')
ylabel('|H_{21}(w)|')
title 'Direct Receptance H_{21}'
grid

subplot(2,2,4)
plot(w/(2*pi), H22)
xlabel('Frequency (Hz)')
ylabel('|H_{22}(w)|')
title 'Direct Receptance H_{22}'
grid

%% Modal receptance matrices

w=linspace(0, 2*pi*200, 1000);
H=zeros(N,N,length(w));

K_mod= Phi.'*K*Phi;
C_mod= Phi.'*C*Phi;
C_mod_diag=diag(diag(C_mod)); 
M_mod= Phi.'*M*Phi;

for n=1:length(w)
    Z_mod=(-w(n)^2*M_mod+1i*w(n)*C_mod_diag+K_mod);
    H(:,:,n)=inv(Z_mod);    
end

for j=1:length(w)
    
H11(j)=abs(H(1,1,j));
H12(j)=abs(H(1,2,j));
H21(j)=abs(H(2,1,j));
H22(j)=abs(H(2,2,j));

end

figure

semilogy(w/(2*pi), H11)
hold on

%semilogy(w/(2*pi), H12)     %Via M_Meth these are 0 and thus are not
%represented


%semilogy(w/(2*pi), H21)


semilogy(w/(2*pi), H22)
xlabel('Frequency (Hz)')
ylabel('|H_{mod}_{j}(w)|')
title 'Modal Receptance H_{11} and H_{22}'
legend ('H11','H22')

ylim ([1e-7 1e-4])
grid

%% Excitació

F=[300 -200].';
we=2*pi*50;

F_mod=Phi.'*F;
K_mod= Phi.'*K*Phi;
C_mod= Phi.'*C*Phi;
C_mod_diag=diag(diag(C_mod));
M_mod= Phi.'*M*Phi;

      
    Z_mod=(-we^2*M_mod+1i*we*C_mod_diag+K_mod);
    Q_mod=Z_mod\F_mod;
    

abs(F_mod)
H=abs(inv(Z_mod))
Q_modul=abs(Q_mod)
Q_angle=rad2deg(angle(Q_mod));
    

toc
