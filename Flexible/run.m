clc
clear
close all

addpath('GFunction');

global  m1 m2 m3 l1 l2 l3 l4 g A_1 A_2 A_3 E1 E2 E3 I1 I2 I3 Rho1 Rho2 Rho3

% l1=0.5;
% l2=1;
% l3=0.75;
% l4=0.45;
% % m1=0.6107;
% % m2=1.2214;
% % m3=0.9161;
% g=0;
% Rho1=2.7143E3;
% Rho2=Rho1;
% Rho3=Rho1;
% E1=71.7e9;
% E2=E1; E3=E1;
% A_1=0.75*3*10^-4;
% A_2=A_1; A_3=A_1;
% I1=0.1055*10^-8;
% I2=I1; I3=I1;


% A_1=7.854E-5;
% A_2=A_1;
% Rho1=2.77E3;
% Rho2=Rho1;
% E1=1e9;
% E2=0.5E8;
% I1=4.909E-10;
% I2=I1;
% m1=A_1*Rho1*l1;
% m2=A_2*Rho2*l2;

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



dt=0.01;
timeSpan=0:dt:1.4;
tsize=length(timeSpan);


% q_0=zeros(12,1);
q_0=[[pi/4;-0.6286 ;-0.3179];zeros(9,1)];
dq_0=zeros(12,1);
mu_0=0;

options=odeset('maxstep',1e-4);
%%
t0=clock;
[~,zGM]=ode45(@GM_dynamics,timeSpan,[q_0;dq_0],options);
t1=clock;
timeGM=etime(t1,t0);
disp(['GM sim time: ' num2str(timeGM) '(s)'])


%%

%
GM.q1=zGM(:,1);
GM.q2=zGM(:,2);
GM.q3=zGM(:,3);

GM.q4=zGM(:,4);
GM.q5=zGM(:,5);
GM.q6=zGM(:,6);
GM.q7=zGM(:,7);
GM.q8=zGM(:,8);
GM.q9=zGM(:,9);
GM.q10=zGM(:,10);
GM.q11=zGM(:,11);
GM.q12=zGM(:,12);

GM.dq1=zGM(:,13);
GM.dq2=zGM(:,14);
GM.dq3=zGM(:,15);

GM.dq4=zGM(:,16);
GM.dq5=zGM(:,17);
GM.dq6=zGM(:,18);
GM.dq7=zGM(:,19);
GM.dq8=zGM(:,20);
GM.dq9=zGM(:,21);
GM.dq10=zGM(:,22);
GM.dq11=zGM(:,23);
GM.dq12=zGM(:,24);
%% 
% %constraint
AE_GM = AEfunc(l1,l2,l3,l4,GM.q1,GM.q2,GM.q3,GM.q5,GM.q8,GM.q11);

% % Energy
% E_GM=Efunc(GM.dq1,GM.dq2,GM.dq3,g,l1,l2,l3,m1,m2,m3,GM.q1,GM.q2,GM.q3);
% 
% 
% ER_GM=(E_GM-E_GM(1))/E_GM(1);
% 
% %%
% figure
% hold on; grid on
% % plot(timeSpan,ER_IM,'r-')
% % plot(timeSpan,ER_AM,'g-')
% % plot(timeSpan,ER_EM,'b-')
% plot(timeSpan,E_GM,'k-')
% legend('GM')
% xlabel('time(s)')
% title('Mechanical Energy Error')
% 
figure
hold on; grid on
% plot(timeSpan,AE_IM,'r-','linewidth',1.5)
% plot(timeSpan,AE_AM,'g.','linewidth',2.5)
% plot(timeSpan,AE_EM,'b-.','linewidth',1.5)
plot(timeSpan,AE_GM,'k--','linewidth',1.5)
legend('GM')
xlabel('time(s)')
grid minor
title('Constraint Error')
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
% 
% figure
% % A=openfig('untitled.fig');
% hold on
% sl=(l1)*cos(GM.q1)+(l2)*cos(GM.q2);
% plot(timeSpan,sl,'-k','LineWidth',1.5);
% % legend('Flex','Rigid');
% grid minor
% xlabel('time(s)')
% ylabel('slider position(m)')
% ylim([0.1 0.5])

%%
yp=sin(GM.q1).*(l1 + GM.q5) + sin(GM.q2).*(l2 + GM.q8);
xp=cos(GM.q1).*(l1 + GM.q5) + cos(GM.q2).*(l2 + GM.q8) ;

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
