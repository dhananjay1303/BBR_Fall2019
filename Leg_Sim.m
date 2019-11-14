clear
close all

%% Import Experimental Data
global Xf Zf
load('FootPlacement.mat')
Xf = fit(timex, 0.8*(x/100),'smoothingspline');
Zf = fit(timez, 0.8*(z/100),'smoothingspline');

%% System Parameter Definitions
global lh lk la lmh lmk lma Ih Ik Ia Mh Mk Ma qstar count bigCount

%Lengths of Leg Sections (meters)
lh = 0.424; %Upper leg length               
lk = 0.425; %Lower leg length               
la = 0.034; %Foot length from ankle to cg   

%Center of mass location from proximal joint (meters)
lmh = 0.184;
lmk = 0.184;
lma = 0.034;

%Moment of inertia of the leg segments (kgm^2)
Ih = 0.113;
Ik = 0.051;
Ia = 0.001;

%Masses of leg segments (kilograms)
Mh = 7.000;
Mk = 3.255;
Ma = 1.015;

%Joint Limits for Comfort (radians)
qh_star = deg2rad(00);
qk_star = deg2rad(10);
qa_star = deg2rad(90);
q_star = [qh_star;qk_star;qa_star];

%Initial State Vector
qh_0 = deg2rad(-5);  % 10
qk_0 = deg2rad(10); %-10
qa_0 = deg2rad(75); %110

%Initial End Point Position
x0 = lh*sin(qh_0) + lk*sin(qh_0-qk_0) + la*sin(qh_0-qk_0+qa_0);
z0 = -lh*cos(qh_0) - lk*cos(qh_0-qk_0) - la*cos(qh_0-qk_0+qa_0);

%Comfortable Leg Positions
qstar = deg2rad([0;10;90]);
t = 0:0.001:timex(end);

count = 4;
bigCount = 0;

% Run ODE solver for desired trajectory
[TOUT_END, ENDPOINT] = ode45(@v, [0 timex(end)], [x0;z0]);

%% Minimum Joint Velocity Squared

[TOUT_JOINT_VEL, JOINTSPACE_VEL] = ode45(@state_deriv_vel,[0 timex(end)],[qh_0;qk_0;qa_0]);

figure(1)
plot(ENDPOINT(:,1),ENDPOINT(:,2)) %Plot desired trajectory
drawnow
hold on

for i = 1:length(TOUT_JOINT_VEL)
    plotLeg(JOINTSPACE_VEL(i,:))  %Plot leg position over time
    if i > 1
        pause(TOUT_JOINT_VEL(i) - TOUT_JOINT_VEL(i-1))
    end
end

%% Minimum Kinetic Energy

[TOUT_JOINT_KE, JOINTSPACE_KE] = ode45(@state_deriv_KE,[0 timex(end)],[qh_0;qk_0;qa_0]);

figure(2)
plot(ENDPOINT(:,1),ENDPOINT(:,2)) %Plot desired trajectory
drawnow limitrate
hold on

for i = 1:2:length(TOUT_JOINT_KE)
    plotLeg(JOINTSPACE_KE(i,:))   %Plot leg position over time
    if i > 1
        pause(TOUT_JOINT_KE(i) - TOUT_JOINT_KE(i-1))
    end
end

%% Joint Velcotiy Plot (both cases)

q_dot_vel = zeros(3,length(TOUT_JOINT_VEL));
q_dot_KE = zeros(3,length(TOUT_JOINT_KE));

for i = 1:length(TOUT_JOINT_VEL)
    %Solve for joint velocities in the minimum velcoity case
    q_dot_vel(:,i) = state_deriv_vel(TOUT_JOINT_VEL(i),transpose(JOINTSPACE_VEL(i,:)));
end

for i = 1:length(TOUT_JOINT_KE)
    %Solve for joint velocities in the minimum kinetic energy case
    q_dot_KE(:,i) = state_deriv_KE(TOUT_JOINT_KE(i), transpose(JOINTSPACE_KE(i,:)));
end

q_dot_cost_vel = sum(sum(q_dot_vel.^2))/length(q_dot_vel);
q_dot_cost_KE  = sum(sum(q_dot_KE .^2))/length(q_dot_KE );

%Plot on the same figure for comparison
figure(3)
plot(TOUT_JOINT_VEL,q_dot_vel(1,:).^2 + q_dot_vel(2,:).^2 + q_dot_vel(3,:).^2)
hold on
plot(TOUT_JOINT_KE, q_dot_KE(1,:).^2 + q_dot_KE(2,:).^2 + q_dot_KE(3,:).^2)
max_time = max(max(TOUT_JOINT_VEL),max(TOUT_JOINT_KE));
xlim([0 max_time])
plot([0 max_time],[0 0],'k')
legend('Vel','KE')
ylabel('Joint Velocity Squared [(°/s)^2]')
xlabel('Time')
grid on

%% Kinetic Energy Plot (both cases)

kinetic_vel = zeros(1,length(TOUT_JOINT_VEL));
kinetic_KE = zeros(1,length(TOUT_JOINT_KE));

for i = 1:length(TOUT_JOINT_VEL)
    %Solve for kinetic energy in the minimum joint velocity case
    q_dot = state_deriv_vel(TOUT_JOINT_VEL(i),transpose(JOINTSPACE_VEL(i,:)));
    kinetic_vel(i) = (1/2)*transpose(q_dot)*H(JOINTSPACE_VEL(i,:))*q_dot;
end

for i = 1:length(TOUT_JOINT_KE)
    %Solve for kinetic energy in the minimum kinetic energy case
    q_dot = state_deriv_KE(TOUT_JOINT_KE(i), transpose(JOINTSPACE_KE(i,:)));
    kinetic_KE(i) = (1/2)*transpose(q_dot)*H(JOINTSPACE_KE(i,:))*q_dot;
end

%Plot on the same figure for comparison
figure(4)
plot(TOUT_JOINT_VEL,kinetic_vel,'r--')
hold on
plot(TOUT_JOINT_KE, kinetic_KE, 'r')
max_time = max(max(TOUT_JOINT_VEL),max(TOUT_JOINT_KE));
xlim([0 max_time])
plot([0 max_time],[0 0],'k')
grid on
legend('Min Velocity','Min KE')

KE_cost_vel = sum(kinetic_vel)/length(kinetic_vel);
KE_cost_KE  = sum(kinetic_KE )/length(kinetic_KE );

fprintf('\n\t\t\t\t\tMinimum Velocity\tMinimum Kinetic Energy\n')
fprintf('\tVelocity Cost\t\t %.2f \t\t\t\t   %.2f\n',q_dot_cost_vel,q_dot_cost_KE)
fprintf('\tKE Cost      \t\t %.2f \t\t\t\t   %.2f\n\n',KE_cost_vel,KE_cost_KE)

%% Functions

function J = Jacobian(q)
    qh = q(1); qk = q(2); qa = q(3);

    global lh lk la
    
    J = zeros(2,3);
    J(1,1) = lh*cos(qh) + lk*cos(qh-qk) + la*cos(qh-qk+qa);
    J(1,2) = -lk*cos(qh-qk) - la*cos(qh-qk+qa);
    J(1,3) = la*cos(qh-qk+qa);
    J(2,1) = lh*sin(qh) + lh*sin(qh-qk) + la*sin(qh-qk+qa);
    J(2,2) = -lh*sin(qh-qk) - la*sin(qh-qk+qa);
    J(2,3) = la*sin(qh-qk+qa);
end

function massMat = H(q)
    qh = q(1); qk = q(2); qa = q(3);
    
    global lmh lmk lma Ih Ik Ia lh lk Mh Mk Ma
    
    massMat = zeros(3,3);
    
    massMat(1,1) = Ih + Mh*lmh^2 + Ik + Mk*(lh^2 + lmk^2 + 2*lh*lmk*cos(qk))...
                 + Ia + Ma*(lh^2 + lk^2 + lma^2 + 2*lh*lk*cos(qk) + 2*lk*lma*cos(qa) + 2*lh*lma*cos(qk + qa));
    massMat(1,2) = Ik + Mk*(lmk^2 + lh*lmk*cos(qk)) + Ia...
                 + Ma*(lk^2 + lma^2 + lh*lk*cos(qk) + 2*lk*lma*cos(qa) + lh*lma*cos(qk + qa));
    massMat(1,3) = Ia + Ma*(lma^2 + lk*lma*cos(qa)) + lh*lma*cos(qk + qa);
    massMat(2,1) = massMat(1,2);
    massMat(2,2) = Ik + Mk*lmk^2 + Ia + Ma*(lk^2 + lma^2 + 2*lk*lma*cos(qa));
    massMat(2,3) = Ia + Ma*(lma^2 + lk*lma*cos(qa));
    massMat(3,1) = massMat(1,3);
    massMat(3,2) = massMat(2,3);
    massMat(3,3) = Ia + Ma*lma^2;
end

%-- Experimental Velocity Data --%
function velocity = v(t,~)
    global Xf Zf
    
    velocity = [differentiate(Xf,t); differentiate(Zf,t)];
end

%-- State Derivative: Minimum Joint Velocity --%
function q_dot = state_deriv_vel(t,q)
    
    global qstar Xf Zf

    xz_dot = [differentiate(Xf,t);differentiate(Zf,t)];
    
    %Calculate pseudo inverse of J to find angular rates
    J = Jacobian(q);
    JT = transpose(J);
    J_pseudo = JT/(J*JT);
    
    q_dot_1 = J_pseudo*xz_dot;
    
    %-- Optimize to avoid uncomfortable joint positions --%
    X = 0; %Rate of convergence
    
    %Project solution onto null space
    n_null = null(J);
    q_dot_null = -X*dot((q-qstar)/norm(q-qstar), n_null)*n_null;
    
    %Determine the final q_dot solution
    q_dot = q_dot_1 + q_dot_null;
    
end

%-- State Derivative: Minimum Kinetic Energy --%
function q_dot = state_deriv_KE(t,q)
    
    global qstar Xf Zf

    xz_dot = [differentiate(Xf,t);differentiate(Zf,t)];
    
    %Calculate pseudo inverse of J to find angular rates
    J = Jacobian(q);
    JT = transpose(J);
    M = H(q);
    M_inv = inv(M);
    
    J_pseudo = M_inv*JT*inv((J*M_inv*JT));
    
    q_dot_1 = J_pseudo*xz_dot;
    
    %-- Optimize to avoid uncomfortable joint positions --%
    if t < 0.3
        X = 10; %Rate of convergence
    else
        X = 0; %Only enforcing cost function in sensitive period of movement
    end
    
    %Project solution onto null space
    n_null = null(J);
    q_dot_null = -X*dot((q-[q(1);q(2);qstar(3)])/norm(q-[q(1);q(2);qstar(3)]), n_null)*n_null;
    
    %Determine the final q_dot solution
    q_dot = q_dot_1 + q_dot_null;
    
end

%-- Plotting Arm Position --%
function plotLeg(q)

    global lh lk la
    qh = q(1); qk = q(2); qa = q(3);

    %shoulder position
    xh = 0;
    yh = 0;

    %knee position
    xk = xh + lh*sin(qh);
    yk = yh - lh*cos(qh);

    %ankle position
    xa = xk + lk*sin(qh-qk);
    ya = yk - lk*cos(qh-qk);

    %end point
    x = xa + la*sin(qh-qk+qa);
    y = ya - la*cos(qh-qk+qa);

    plot([xh xk xa x],[yh yk ya y],'color',0*[1 1 1])
    drawnow
end