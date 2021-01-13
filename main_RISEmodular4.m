%% MAINSCRIPT RISE-based modular adaptive controller with update law thetahatdot=gamma*Yd'*tanh(e);
clc; clear; close all;

global Ks gamma alpha1 alpha2 beta E0 


%%%%%%%%%%%%%%%%%%%Tunable Controller Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ks=9; %Constant in mu
beta=20;%Constant multiplying sgn(E) in the expression for nudot
alpha1=1; %Constant in the definition of the error signal 'E'
alpha2=1; %Constant in the definition of the error signal 'r'
gamma=diag([300; 2; 20; 20 ;20]);%Learning rate matrix in the adaptive update law

%% Augmented Initial State:
% Initial State:
q10=3; q20=4.5; q1dot0=0; q2dot0=0;
x0=[q10;q20;q1dot0;q2dot0];
%Calculating the initial value of error signal E for use in the control law
t=0; q0=[q10;q20]; qdot0=[q1dot0;q2dot0];
qd10=cos(0.5*t); qd20=2*cos(t); qd1dot0=-0.5*sin(0.5*t); qd2dot0=-2*sin(t);
qd0=[qd10;qd20];  qd_dot0=[qd1dot0;qd2dot0];
e0=qd0-q0; 
edot0=qd_dot0-qdot0;
E0= edot0 + alpha1*e0;

% Initial Values of the adaptive estimates: 
%(We pick these based on our expectations of the actual parameter values)
p1hat0=8; p2hat0=0.7; p3hat0=0.5; fd1hat0=8; fd2hat0=2; 
thetahat0=[p1hat0;p2hat0;p3hat0;fd1hat0;fd2hat0];

%Initial value of nu appearing in the control law which is to be solved for by numerical integration
nu0=[0;0];

%Augmented Initial State
X0=[x0;thetahat0; nu0];

%% Simulation time
t_sim=20;

%% Obtaining the solution to the ode based on the control law specified in ode_compadaptgrad
tspan=[0 t_sim];
options = odeset('AbsTol',1e-3,'RelTol',1e-3);
[time,X_Sol]=ode113(@ode_RISEmodular4,tspan,X0,options);

%% Extracting the trajectories of the states and the estimates from the ODE solution matrix X_Sol:
q1=X_Sol(:,1);
q2=X_Sol(:,2);
q1dot=X_Sol(:,3);
q2dot=X_Sol(:,4);
p1hat=X_Sol(:,5);
p2hat=X_Sol(:,6);
p3hat=X_Sol(:,7);
fd1hat=X_Sol(:,8);
fd2hat=X_Sol(:,9);
nu_1=X_Sol(:,10);
nu_2=X_Sol(:,11);
%% Computing desired trajectories, errors and controls from the ode solution for the plots:

%Desired trajectories:
qd1=cos(0.5*time); qd2=2*cos(time);
qd1dot=-0.5*sin(0.5*time); qd2dot=-2*sin(time);
qd1doubledot=-0.25*cos(0.5*time); qd2doubledot=-2*cos(time);

%Actual Parameter values (Have not been used in the control law)
p1 = 3.473;
p2 = 0.196;
p3 = 0.242;
fd1 = 5.3;
fd2 = 1.1;

%% Errors
e1=qd1-q1; e2=qd2-q2;
e1dot=qd1dot-q1dot; e2dot=qd2dot-q2dot;
p1error= p1-p1hat; p2error= p2-p2hat; p3error= p3-p3hat; fd1error= fd1-fd1hat; fd2error= fd2-fd2hat;
E1=e1dot+alpha1*e1;
E2=e2dot+alpha1*e2;

%% Controls
ydthetahat_1= fd1hat.*qd1dot + p1hat.*qd1doubledot + p2hat.*qd2doubledot + p3hat.*(2.*qd1doubledot.*cos(qd2) + qd2doubledot.*cos(qd2) - qd1dot.*qd2dot.*sin(qd2) - qd2dot.*sin(qd2).*(qd1dot + qd2dot));
ydthetahat_2=p2hat.*(qd1doubledot + qd2doubledot) + fd2hat.*qd2dot + p3hat.*(sin(qd2).*qd1dot.^2 + qd1doubledot.*cos(qd2));
tau1= ydthetahat_1 + (Ks+1)*E1 - (Ks+1)*E0(1) + nu_1;
tau2= ydthetahat_2 + (Ks+1)*E2 - (Ks+1)*E0(2) + nu_2;

%% Plots
% figure;
subplot(2,2,1);
plot(time,q1,'r',time,q2,'b',time,qd1,'r--',time,qd2,'b--','LineWidth',1.5);
xlabel('Time');
ylabel('States');
title('Desired Vs Actual Trajectories');
xlim([0 time(end)]);
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
leg1 = legend('$q_1$','$q_2$','$q_{d1}$','$q_{d2}$');
set(leg1,'Interpreter','latex');
hold on;

% figure;
subplot(2,2,2);
plot(time,e1,'r',time,e2,'b','LineWidth',1.5);
xlabel('Time');
ylabel('Tracking Errors');
title('Tracking Errors Vs Time');
xlim([0 time(end)]);
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
leg3 = legend('$e_1$','$e_2$');
set(leg3,'Interpreter','latex');
hold on;

% figure;
subplot(2,2,3);
plot(time,p1hat,'r',time,p2hat,'b',time,p3hat,'k',time,fd1hat,'m',time,fd2hat,'g','LineWidth',1.5);
xlabel('Time');
ylabel('Adaptive Estimates');
title('Adaptive Estimates Vs Time');
xlim([0 time(end)]);
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
leg2 = legend('$\widehat{p}_1$','$\widehat{p}_2$','$\widehat{p}_3$','$\widehat{f}_{d1}$','$\widehat{f}_{d2}$');
set(leg2,'Interpreter','latex');
hold on;

% figure;
subplot(2,2,4);
plot(time,p1error,'r',time,p2error,'b',time,p3error,'k',time,fd1error,'m',time,fd2error,'g','LineWidth',1.5);
xlabel('Time');
ylabel('Actual value - Adaptive Estimate');
title('Parameter Estimate Errors');
xlim([0 time(end)]);
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
leg2 = legend('$p_1-\widehat{p}_1$','$p_2-\widehat{p}_2$','$p_3-\widehat{p}_3$','$f_{d1}-\widehat{f}_{d1}$','$f_{d2}-\widehat{f}_{d2}$');
set(leg2,'Interpreter','latex');
hold on;


figure;
plot(time,tau1,'r',time,tau2,'b','LineWidth',1.5);
xlabel('Time');
ylabel('Controls');
title('Controls Vs Time');
xlim([0 time(end)]);
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.3;
ax.FontSize = 16;
ax.LineWidth = 1.4;
leg4 = legend('$\tau_1$','$\tau_2$');
set(leg4,'Interpreter','latex');
hold on;

%% The following lines are useful for tuning the gains:
maxlink1torque=max(abs(tau1))
maxlink2torque=max(abs(tau2))