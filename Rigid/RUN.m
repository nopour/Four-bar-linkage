clc
clear
close all

addpath('GFunction');

global  l1 l2 l3 l4 m1 m2 m3 g Rho1 Rho2 A_1 A_2 E1 E2 I1 I2

% l1=1;
% l2=2;
% l3=1.5;
% l4=0.9;
% m1=0.6107;
% m2=1.2214;
% m3=0.9161;
% g=0;


l1=1;
l2=2;
l3=1.5;
l4=0.9;
g=-9.8;
Rho1=2.7143E3;
Rho2=Rho1;
Rho3=Rho1;
E1=71.7e9;
% E1=3e9;
E2=E1; E3=E1;
A_1=0.75*3*10^-4;
A_2=A_1; A_3=A_1;
I1=0.1055*10^-8;
I2=I1; I3=I1;

m1=A_1*Rho1*l1;
m2=A_2*Rho2*l2;
m3=A_3*Rho3*l3;


% l1=5;
% l2=11;
% l3=10.5;
% l4=10;
% g=-32;
% Rho1=0.000725;
% Rho2=Rho1;
% Rho3=Rho1;
% E1=30E6;
% E2=E1; E3=E1;
% A_1=0.25;
% A_2=A_1; A_3=A_1;
% I1=0.0208;
% I2=I1; I3=I1;
% 
% m1=A_1*Rho1*l1;
% m2=A_2*Rho2*l2;
% m3=A_3*Rho3*l3;


dt=0.01;
timeSpan=0:dt:1.4;
tsize=length(timeSpan);


% q_0=[0;1.23978 ;0.1361363267289];

q_0=[pi/4;-0.6286 ;-0.3179];
% q_0=[pi/4;-1.9804 ;-2.2910];

% q_0=[0;0;0;0;0;0;0;0];
dq_0=[0;0;0];
mu_0=0;

options=odeset('maxstep',1e-4);
%%
t0=clock;
[~,zGM]=ode45(@AM_dynamics,timeSpan,[q_0;dq_0],options);
t1=clock;
GMT=etime(t1,t0);
disp(['GM sim time: ' num2str(GMT) '(s)'])


%%

%
GM.q1=zGM(:,1);
GM.q2=zGM(:,2);
GM.q3=zGM(:,3);
GM.dq1=zGM(:,4);
GM.dq2=zGM(:,5);
GM.dq3=zGM(:,3);
%% 
%constraint
AE_GM = AEfunc(l1,l2,l3,l4,GM.q1,GM.q2,GM.q3);


figure
hold on; grid on
% plot(timeSpan,AE_IM,'r-','linewidth',1.5)
% plot(timeSpan,AE_AM,'g.','linewidth',2.5)
% plot(timeSpan,AE_EM,'b-.','linewidth',1.5)
plot(timeSpan,AE_GM,'color','#77AC30','linewidth',2)
legend('Agumented Method')
xlabel('time(s)')
grid minor
title('$ Constraint $ $ Error $','Interpreter','latex')
% ylabel('$ Y [m]$','Interpreter','latex')
xlabel('$ Time $ $ [s]$','Interpreter','latex')
%%
figure
hold on; grid on
% plot(timeSpan,(IM.q1)*180/pi,'r-','linewidth',2.5)
% plot(timeSpan,(AM.q1)*180/pi,'g--','linewidth',2)
% plot(timeSpan,(EM.q1)*180/pi,'b-','linewidth',1.5)
plot(timeSpan,(GM.q1)*180/pi,'k-','linewidth',1)
legend('GM')
xlabel('time(s)')
title('\theta_1')

%%
figure
hold on; grid on
% plot(timeSpan,(IM.q2)*180/pi,'r-','linewidth',2.5)
% plot(timeSpan,(AM.q2)*180/pi,'g--','linewidth',2)
% plot(timeSpan,(EM.q2)*180/pi,'b-','linewidth',1.5)
plot(timeSpan,(GM.q2)*180/pi,'k-','linewidth',1)
legend('GM')
xlabel('time(s)')
title('\theta_2')
ylabel('Angle (deg)')

%%

yp=l1*sin(GM.q1)+l2*sin(GM.q2);
xp=l1*cos(GM.q1)+l2*cos(GM.q2);
figure
plot(timeSpan,yp,'k-','linewidth',1.5)
grid minor
ylabel('$ y [m]$','Interpreter','latex')
xlabel('$ Time [s]$','Interpreter','latex')

figure
plot(timeSpan,xp,'b-','linewidth',1.5)
grid minor
ylabel('$ X [m]$','Interpreter','latex')
xlabel('$ Time [s]$','Interpreter','latex')
