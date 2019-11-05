clear
close all

%% Import Experimental Data
load('FootPlacement.mat')
Xf = fit(time_x-time_x(1), X,'poly8');
Zf = fit(time_z-time_z(1), Z,'poly4');

%Create ppolynomials
global PX_dot PZ_dot PX_double_dot PZ_double_dot PX_triple_dot PZ_triple_dot
PX = [Xf.p1 Xf.p2 Xf.p3 Xf.p4 Xf.p5 Xf.p6 Xf.p7 Xf.p8 Xf.p9];
PZ = [Zf.p1 Zf.p2 Zf.p3 Zf.p4 Zf.p5];
PX_dot = polyder(PX);
PZ_dot = polyder(PZ);
PX_double_dot = polyder(PX_dot);
PZ_double_dot = polyder(PZ_dot);
PX_triple_dot = polyder(PX_double_dot);
PZ_triple_dot = polyder(PZ_double_dot);

%% System Parameter Definitions
global lh lk la qstar q_dot_vec t_vec

q_dot_vec = [];

%Lengths of Leg Sections (not set yet)
lh = 0.6; %Upper leg length               (m)
la = 0.6; %Lower leg length               (m)
lk = 0.1; %Foot length from ankle to cg   (m)

%Joint Limits for Comfort
qh_star = deg2rad(00);
qk_star = deg2rad(30);
qa_star = deg2rad(45);
q_star = [qh_star;qk_star;qa_star];

%Initial State Vector
qh_0 = deg2rad( 10);
qk_0 = deg2rad(-10);
qa_0 = deg2rad(110);

%Comfortable Leg Positions
qstar = deg2rad([0;10;90]);

%% Functions

function J = Jacobian(q)
    qh = q(1); qk = q(2); qa = q(3);

    global lh lk la
    
    J = zeros(2,3);
    J(1,1) = lh*cos(qh) + lk*cos(qh-qk) + la*cos(qh-qk+qa);
    J(1,2) = -lk*cos(qh-qk) - la*cos(qh-qk+qa);
    J(1,3) = la*cos(qh-qk+qa);
    J(2,1) = -lh*sin(qh) - lh*sin(qh-qk) - la*sin(qh-qk+qa);
    J(2,2) = lh*sin(qh-qk) + la*sin(qh-qk+qa);
    J(2,3) = -la*sin(qh-qk+qa);
end

function q_dot = state_deriv(t,q)
    
    global PX_dot PZ_dot PX_double_dot PZ_double_dot PX_triple_dot PZ_triple_dot q_dot_vec
    
    if length(q_dot_vec) < 4
        jerk = [0;0;0]; %not enough points to approximate jerk yet
    else
        T4 = t_vec(end-3:end); %extract last 4 time points
        
        %-- Fit last 4 q_dot points to a polynomial --%
        q_dot_hip = fit(T4,q_dot_vec(1,end-3:end),'poly3');
        q_dot_knee = fit(T4,q_dot_vec(2,end-3:end),'poly3');
        q_dot_ankle = fit(T4,q_dot_vec(3,end-3:end),'poly3');
        
        %-- Find the second derivatives to get jerk --%
        [~, jerk_hip] = differentiate(q_dot_hip,T4);
        [~, jerk_knee] = differentiate(q_dot_knee,T4);
        [~, jerk_ankle] = differentiate(q_dot_ankle,T4);
        
        jerk = [jerk_hip;jerk_knee;jerk_ankle];
    end
    
    %Determine velocity of end point
    xz_dot = [polyval(PX_dot,t); polyval(PZ_dot,t)];
    xz_double_dot = [polyval(PX_double_dot,t); polyval(PZ_double_dot,t)];
    xz_triple_dot = [polyval(PX_triple_dot,t); polyval(PZ_triple_dot,t)];

    %Calculate pseudo inverse of J to find angular rates
    J = Jacobian(q);
    JT = transpose(J);
    JJT_inv = inv(J*JT);
    
    %calculate q triple dot
    
    %-- Optimize to avoid uncomfortable joint positions --%
    X = 1; %Rate of convergence
    
    %Project solution onto null space
    n_null = null(J);
    q_dot_null = -X*dot(jerk/norm(jerk), n_null)*n_null;
    
    %Determine the final q_dot solution
    q_dot = q_dot_1 + q_dot_null;
end