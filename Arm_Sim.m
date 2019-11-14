clear
close all

%% Definition of Physical Parameters
%Subscript s refers to shoulder
%Subscript e refers to elbow
%Subscript h refers to hand

global ms me mh ls le lh lms lme lmh Is Ie Ih qstar count bigCount RUN

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

%Define target "comfortable" position for each joint
qstar = 0.5*[qs_plus+qs_minus;qe_plus+qe_minus;qh_plus+qh_minus];

%% State Definition

%Initial joint angle vector
qs0 = deg2rad(0); %0
qe0 = deg2rad(120); %120
qh0 = deg2rad(75); %75
q0 = [qs0;qe0;qh0];


x0 = ls*cos(qs0) + le*cos(qs0+qe0) + lh*cos(qs0+qe0+qh0);
y0 = ls*sin(qs0) + le*sin(qs0+qe0) + lh*sin(qs0+qe0+qh0);

count = 0;
bigCount = 0;

figure(1)

%Calculate desired trajectory
[~, XY] = ode45(@v, [0 1], [x0;y0]);
plot(XY(:,1),XY(:,2),'b','Linewidth',2)
drawnow
hold on

RUN = 1;
%Calculate minimum joint velocity trajectory
[TOUT_JV, JOINTSPACE_JV] = ode45(@qdot, [0 1], q0); %solve for end point over time

figure(2)
plot(XY(:,1),XY(:,2),'b','Linewidth',2)
drawnow
hold on

RUN = 2; bigCount = 0; count = 0;
%Calculate minimum kinetic energy trajectory
[TOUT_KE, JOINTSPACE_KE] = ode45(@qdot, [0 1], q0); %solve for end point over time

%% Evaluate Performance Indices (Joint Velocity and Kinetic Energy)

n_JV = length(TOUT_JV);
n_KE = length(TOUT_KE);

RUN = 3;
q_dot_cost_JV = 0;
KE_cost_JV = 0;
for i = 1:length(TOUT_JV)
    t = TOUT_JV(i);
    q = JOINTSPACE_JV(i,:)';
    q_dot_JV = qdot(t,q);
    q_dot_cost_JV = q_dot_cost_JV + norm(q_dot_JV)^2/n_JV;
    KE_cost_JV = KE_cost_JV + (1/2)*transpose(q_dot_JV)*massMat(q)*q_dot_JV/n_JV;
end

RUN = 4;
q_dot_cost_KE = 0;
KE_cost_KE = 0;
for i = 1:length(TOUT_KE)
    t = TOUT_KE(i);
    q = JOINTSPACE_KE(i,:)';
    q_dot_KE = qdot(t,q);
    q_dot_cost_KE = q_dot_cost_KE + norm(q_dot_KE)^2/n_KE;
    KE_cost_KE = KE_cost_KE + (1/2)*transpose(q_dot_KE)*massMat(q)*q_dot_KE/n_KE;
end
    
fprintf('\n\t\t\t\t\tMinimum Velocity\tMinimum Kinetic Energy\n')
fprintf('\tVelocity Cost\t\t %.2f \t\t\t\t   %.2f\n',q_dot_cost_JV,q_dot_cost_KE)
fprintf('\tKE Cost      \t\t %.2f \t\t\t\t   %.2f\n',KE_cost_JV,KE_cost_KE)

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

%% Cost Function Derivative

function delV = grad_V(q)

global qstar
x = q(1); y = q(2); z = q(3);
delV = 2*[x-qstar(1);y-qstar(2);z-qstar(3)];

end

%% State Propagation Function

function q_dot_star = qdot(t,q)

    global count bigCount RUN
    T = 1; %duration of movement
    A = 0.55; %amplitude of movement

    %Bell shaped velocity profile
    velocity = ((t/T)^2 - 2*(t/T) + 1)*30*A*(t/T)^2/T;
    xydot = velocity*[cosd(20);sind(20)];
    
    %Minimization parameters for kinetic energy
    JT = transpose(J(q));
    M = massMat(q);
    M_inv = inv(M);

    if RUN == 1 || RUN == 3
        J_pseudo = JT/(J(q)*JT);                   % minimize velocity
    else
        J_pseudo = M_inv*JT*inv((J(q)*M_inv*JT)); % minimize kinetic energy
    end
    
    q_dot_star = J_pseudo*xydot;

% Steepest descent optimization to avoid uncomfortable positions
    X = 1.8; %arbitrary, determines the rate of convergence 
    n = -grad_V(q)/norm(grad_V(q));

    %Projection of n onto the null space of J
    NS = null(J(q));
    n_null = (dot(n,NS))*NS/norm(NS);

    %Optimization contribution of joint space
    delta_q_dot = -X*dot(grad_V(q),n_null)*n_null;

    %Combining minimization term with optimization term
    q_dot_star = q_dot_star + delta_q_dot;

    %Only plot every 4th position
    if RUN == 1 || RUN == 2
        if count ==4
            plotArm(q)
            count = 0;
        else
            count = count +1;
        end

        bigCount = bigCount + 1;
        if bigCount > 100
            %assert(false)
        end
    end
end

%% Velocity Profile Function

function xydot = v(t,~)
    T = 1;
    A = 0.55;

    %Bell shaped velocity profile
    velocity = ((t/T)^2 - 2*(t/T) + 1)*30*A*(t/T)^2/T;
    
    %Arbitrarily chosen direction
    xydot = velocity*[cosd(20);sind(20)];
end

%% Plotting Arm Position

function plotArm(q)

global ls le lh bigCount
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

plot([xs xe xh x],[ys ye yh y],'color',[1 1 1].*max(0.7-bigCount/100,0))
end