clc, clear all, close all;
tic
mr=216; % [kg]
ms=62784; % [kg]
J=895124; % [kgm^2]

% kf=10e6; % [N/m] 
%cf=1e3; % [Ns/m] 

% kf=28e6; % [N/m] 
% cf=2.8e3; % [Ns/m] 
% 
% kf=90e6; % [N/m] 
% cf=9e3; % [Ns/m] 

ks=1.5e8; % [N/m]
d1=1.646; % [m]
d2=5.45; % [m]

M=diag([mr ms J]);

%% 1st values
%from 3 different spring + dampener combinations to reduce vibrations in
%L9, BCN subway.

kf=10e6; % [N/m] 
cf=1e3; % [Ns/m] 

K=[kf -kf kf*d1;...
   -kf 2*ks+kf -kf*d1;...
   kf*d1 -kf*d1 kf*d1^2+2*ks*d2^2];

C=[cf -cf cf*d1;...
   -cf cf -cf*d1;...
   cf*d1 -cf*d1 cf*d1^2];

N=length(M); % Number of DOFs of the problem.
% Normal modes
%A=M\K;
DOF=length(M);
[Phi, Lambda]=eig(K,M);
fn=sqrt(Lambda)/(2*pi)
%Phi
for i=1:DOF
    Phi_norm(:,i)=Phi(:,i)/Phi(1,i); 
end
Phi_norm
Phi*10^3;

%Receptance

w=linspace(0, 2*pi*200, 1000);
H=zeros(N,N,length(w));

for n=1:length(w)
    Z=(-w(n)^2*M+1i*w(n)*C+K);
    H(:,:,n)=inv(Z);    
end

for j=1:length(w)
    
H21(j)=abs(H(2,1,j));

end

semilogy(w/(2*pi), H21)
hold on





        % w=linspace(0, 2*pi*200, 1000);
        % F=[1 0 0].';
        % Q_dir=zeros(length(w),1);
        % 
        % for n=1:length(w)
        %   for j=1:3
        %       
        %     Z=(-w(n)^2*M+1i*w(n)*C+K);
        %     Q=Z\F;
        %     Q_dir(j,n)=abs(Q(j));
        % 
        %   end
        % end
        % semilogy(w/(2*pi),Q_dir(2,:))%modeshape 2, demana la transfe entre la força i el moviment vertical de la llosa
        % hold on
        % grid

        
fprintf '///////////////////////////////////////'
fprintf 'RIGIDESA 1'

%% 2nd values

kf=28e6; % [N/m] 
cf=2.8e3; % [Ns/m] 

K=[kf -kf kf*d1;...
   -kf 2*ks+kf -kf*d1;...
   kf*d1 -kf*d1 kf*d1^2+2*ks*d2^2];

C=[cf -cf cf*d1;...
   -cf cf -cf*d1;...
   cf*d1 -cf*d1 cf*d1^2];

N=length(M); % Number of DOFs of the problem.
% Normal modes
%A=M\K;
DOF=length(M);
[Phi, Lambda]=eig(K,M);
fn=sqrt(Lambda)/(2*pi)
%Phi
for i=1:DOF
    Phi_norm(:,i)=Phi(:,i)/Phi(1,i); 
end
Phi_norm
Phi*10^3;

%Receptance


w=linspace(0, 2*pi*200, 1000);
H=zeros(N,N,length(w));

for n=1:length(w)
    Z=(-w(n)^2*M+1i*w(n)*C+K);
    H(:,:,n)=inv(Z);    
end

for j=1:length(w)
    
H21(j)=abs(H(2,1,j));

end

semilogy(w/(2*pi), H21)

        % w=linspace(0, 2*pi*200, 1000);
        % F=[1 0 0].';
        % Q_dir=zeros(length(w),1);
        % 
        % for n=1:length(w)
        %   for j=1:3
        %       
        %     Z=(-w(n)^2*M+1i*w(n)*C+K);
        %     Q=Z\F;
        %     Q_dir(j,n)=abs(Q(j));
        % 
        %   end
        % end
        % semilogy(w/(2*pi),Q_dir(2,:))%modeshape 2, demana la transfe entre la força i el moviment vertical de la llosa
        % hold on

fprintf '///////////////////////////////////////'
fprintf 'RIGIDESA 2'
%% 3rd values

kf=90e6; % [N/m] 
cf=9e3; % [Ns/m] 

K=[kf -kf kf*d1;...
   -kf 2*ks+kf -kf*d1;...
   kf*d1 -kf*d1 kf*d1^2+2*ks*d2^2];

C=[cf -cf cf*d1;...
   -cf cf -cf*d1;...
   cf*d1 -cf*d1 cf*d1^2];

N=length(M); % Number of DOFs of the problem.
% Normal modes
%A=M\K;
DOF=length(M);
[Phi, Lambda]=eig(K,M);
fn=sqrt(Lambda)/(2*pi)
%Phi
for i=1:DOF
    Phi_norm(:,i)=Phi(:,i)/Phi(1,i); 
end
Phi_norm
Phi*10^3;

%% Receptance
%from every modeshape 


w=linspace(0, 2*pi*200, 1000);
H=zeros(N,N,length(w));

for n=1:length(w)
    Z=(-w(n)^2*M+1i*w(n)*C+K);
    H(:,:,n)=inv(Z);    
end

for j=1:length(w)
    
H21(j)=abs(H(2,1,j));

end

semilogy(w/(2*pi), H21)

      
      

fprintf '///////////////////////////////////////'
fprintf 'RIGIDESA 3'


toc
