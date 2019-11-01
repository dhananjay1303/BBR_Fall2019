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
global lh lk la qstar

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

function grad_V = cost_deriv(q)

    global qstar
    x = q(1); y = q(2); z = q(3);
    grad_V = 2*[x-qstar(1);y-qstar(2);z-qstar(3)];
    
end

function q_dot = state_deriv(t,q)
    
    global PX_dot PZ_dot PX_double_dot PZ_double_dot PX_triple_dot PZ_triple_dot
    
    %Determine velocity of end point
    xz_dot = [polyval(PX_dot,t); polyval(PZ_dot,t)];
    xz_double_dot = [polyval(PX_double_dot,t); polyval(PZ_double_dot,t)];
    xz_triple_dot = [polyval(PX_triple_dot,t); polyval(PZ_triple_dot,t)];

    %Calculate pseudo inverse of J to find angular rates
    J = Jacobian(q);
    JT = transpose(J);
    JJT_inv = inv(J*JT);
    JTJ_inv = inv(JT*J);
    
    %calculate q triple dot
    
    %Optimize to avoid uncomfortable joint positions
    X = 1; %Rate of convergence
    n = -cost_deriv(q)/norm(cost_deriv(q));
    
    %Project solution onto null space
    NS = null(J);
    n_null = (dot(n,NS))*NS/norm(NS);
    
    %Scale null space projection
    q_dot_null = -X*dot(cost_deriv(q),n_null)*n_null;
    
    %Determine the final q_dot solution
    q_dot = q_dot_1 + q_dot_null;
end