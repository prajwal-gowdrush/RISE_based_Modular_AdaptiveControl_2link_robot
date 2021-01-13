function Xdot=ode_RISEmodular4(t,X)
global Ks gamma alpha1 alpha2 beta E0 

%Assigning states from the input arguments to other variables
q1=X(1); q2=X(2); q1dot=X(3); q2dot=X(4); 
p1hat=X(5); p2hat=X(6); p3hat=X(7); fd1hat=X(8); fd2hat=X(9);
nu=X(10:11,1);
q=[q1;q2]; qdot=[q1dot;q2dot]; 
thetahat=[p1hat;p2hat;p3hat;fd1hat;fd2hat];

%Defining the desired trajectory
qd1=cos(0.5*t); qd2=2*cos(t);
qd1dot=-0.5*sin(0.5*t); qd2dot=-2*sin(t);
qd1doubledot=-0.25*cos(0.5*t); qd2doubledot=-2*cos(t);
qd=[qd1;qd2]; 
qd_dot=[qd1dot;qd2dot];

%Error definitions
e=qd-q; 
edot=qd_dot-qdot;
E= edot + alpha1*e;

%Regression Matrices
Yd =[ qd1doubledot, qd2doubledot, 2*qd1doubledot*cos(qd2) + qd2doubledot*cos(qd2) - qd1dot*qd2dot*sin(qd2) - qd2dot*sin(qd2)*(qd1dot + qd2dot), qd1dot, 0; 0, qd1doubledot + qd2doubledot, sin(qd2)*qd1dot^2 + qd1doubledot*cos(qd2), 0, qd2dot];  

%A term used in the control law
mu= (Ks+1)*E - (Ks+1)*E0 + nu;

%Control Law
tau= Yd*thetahat + mu;

%Actual Parameter values (Have not been used in the control law)
p1 = 3.473;
p2 = 0.196;
p3 = 0.242;
fd1 = 5.3;
fd2 = 1.1;

%Adaptive Update Law
thetahatdot=gamma*Yd'*tanh(e);

%Matrices appearing in the dynamics
M= [ p1+2*p3*cos(q2), p2+p3*cos(q2); p2+p3*cos(q2),  p2]; %Inertia Matrix
V= [ -p3*sin(q2)*q2dot, -p3*sin(q2)*(q1dot+q2dot); p3*sin(q2)*q1dot , 0]; %Centrifugal Coriolis Matrix
Fd=[fd1 , 0;  0 , fd2];

%% Dynamics
Xdot(1,1)= q1dot;
Xdot(2,1)= q2dot;
Xdot(3:4,1)= M\(tau - V*qdot - Fd*qdot); % M\ is inverse of M
Xdot(5:9,1)=thetahatdot;
Xdot(10:11,1)=(Ks+1)*alpha2*E + beta*sign(E);