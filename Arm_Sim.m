clear
close all

%% Definition of Physical Parameters
%Subscript s refers to shoulder
%Subscript e refers to elbow
%Subscript h refers to hand

global ms me mh ls le lh lms lme lmh Is Ie Ih qstar count

%masses of arm segments (kg)
ms = 1.93;
me = 1.52;
mh = 0.52;

%lengths of arm segments (m)
ls = 0.31;
le = 0.34;
lh = 0.18;

%center of mass positions from proximal joint (m)
lms = 0.165;
lme = 0.19;
lmh = 0.055;

%moments of inertia of arm segments (kgm^2)
Is = 0.0141;
Ie = 0.0188;
Ih = 0.0003;

%joint limits
qs_plus = deg2rad(30);
qs_minus = deg2rad(-60);
qe_plus = deg2rad(150);
qe_minus = deg2rad(0);
qh_plus = deg2rad(90);
qh_minus = deg2rad(-75);
qstar = 0.5*[qs_plus+qs_minus;qe_plus+qe_minus;qh_plus+qh_minus];

%% State Definition

%Initial joint angle vector
q0 = deg2rad([0;120;75]);
qs0 = q0(1); qe0 = q0(2); qh0 = q0(3);

count = 0;
figure(1)
hold on
[TOUT, POSITION] = ode45(@qdot, [0 1], q0); %solve for end point over time

%% Inertial Model

function H = massMat(q)

global ms me mh ls le lms lme lmh Is Ie Ih

qe = q(1); qh = q(2);
H = zeros(3,3);

H(1,1) = Is + ms*lms^2 + Ie + me*(ls^2 + lme^2 + 2*ls*lme*cos(qe)) + Ih ...
    + mh*(ls^2 + le^2 +lmh^2 + 2*ls*le*cos(qe) + 2*le*lmh*cos(qh) + 2*ls*lmh*cos(qe+qh));
H(1,2) = Ie + me*(lme^2 + ls*lme*cos(qe)) + Ih + mh*(le^2 + lmh^2 +ls*le*cos(qe) ...
    + 2*le*lmh*cos(qh) + ls*lmh*cos(qe+qh));
H(2,1) = H(1,2);
H(1,3) = Ih + mh*(lmh^2 + 2*le*lmh*cos(qh)) + ls*lmh*cos(qe+qh);
H(3,1) = H(1,3);
H(2,2) = Ie + me*lme^2 + Ih + mh*(le^2 + lmh^2 + 2*le*lmh*cos(qh));
H(2,3) = Ih + mh*(lmh^2 + le*lmh*cos(qh));
H(3,2) = H(2,3);
H(3,3) = Ih + mh*lmh^2;
end

%% Jacobian Model

function Jacobian = J(q)

global ls le lh
qs = q(1); qe = q(2); qh = q(3);

Jacobian = zeros(2,3);
Jacobian(1,1) = -ls*sin(qs) - le*sin(qs+qe) - lh*sin(qs+qe+qh);
Jacobian(1,2) = -le*sin(qs+qe) - lh*sin(qs+qe+qh);
Jacobian(1,3) = -lh*sin(qs+qe+qh);
Jacobian(2,1) =  ls*cos(qs) + le*cos(qs+qe) + lh*cos(qs+qe+qh);
Jacobian(2,2) =  le*cos(qs+qe) + lh*cos(qs+qe+qh);
Jacobian(2,3) = lh*cos(qs+qe+qh);

end

%% Cost Function and its derivative

function delV = grad_V(q)

global qstar
x = q(1); y = q(2); z = q(3);
delV = 2*[x-qstar(1);y-qstar(2);z-qstar(3)];

end

function cost = V(q)

global qstar
x = q(1); y = q(2); z = q(3);
cost = (x-qstar(1))^2 + (y - qstar(2))^2 + (z - qstar(3))^2;

end

%% State Propagation Function

function q_dot_star = qdot(t,q)

global count
T = 1; %duration of movement
A = 0.4; %amplitude of movement

%Bell shaped velocity profile
velocity = ((t/T)^2 - 2*(t/T) + 1)*30*A*(t/T)^2/T;
xydot = velocity*[cosd(20);sind(20)];

%Minimization parameters for kinetic energy
JT = transpose(J(q));
M = massMat(q);
M_inv = inv(M);

J_pseudo = M_inv*JT*inv((J(q)*M_inv*JT));
q_dot_star = J_pseudo*xydot;

%Steepest descent optimization to avoid uncomfortable positions
X = 0.5; %calculated by hand
n = -grad_V(q)/norm(grad_V(q));
delta_q_dot = -X*dot(grad_V(q),n)*n;
q_dot_star = q_dot_star + delta_q_dot;

if count ==4
    plotArm(q)
    count = 0;
    pause
else
    count = count +1;
end

end

%% Plotting Arm Position

function plotArm(q)

global ls le lh
qs = q(1); qe = q(2); qh = q(3);

%shoulder position
xs = 0;
ys = 0;

%elbow position
xe = ls*cos(qs);
ye = ls*sin(qs);

%wrist position
xh = ls*cos(qs) + le*cos(qs+qe);
yh = ls*sin(qs) + le*sin(qs+qe);

%end point
x = ls*cos(qs) + le*cos(qs+qe) + lh*cos(qs+qe+qh);
y = ls*sin(qs) + le*sin(qs+qe) + lh*sin(qs+qe+qh);

plot([xs xe xh x],[ys ye yh y])
end