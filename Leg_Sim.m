clear
close all

%% System Parameter Definitions
global lh lk la qstar

%Lengths of Leg Sections (not set yet)
lh = 1; %Upper leg length               (m)
la = 1; %Lower leg length               (m)
lk = 1; %Foot length from ankle to cg   (m)

%Joint Limits for Comfort
qh_star = deg2rad(00);
qk_star = deg2rad(30);
qa_star = deg2rad(45);
q_star = [qh_star;qk_star;qa_star];

%Initial State Vector
qh_0 = deg2rad( 10);
qk_0 = deg2rad(-10);
qa_0 = deg2rad(110);

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

function xy_dot = cartesian_deriv(t,~)
    %curve-fit
    xy_dot = t;
end

function grad_V = cost_deriv(q)
    global qstar
    %TBD
end

function q_dot = state_deriv(t,q)
    
    T = 1; %Duration of movement
    
    %Determine velocity of end point
    xy_dot = cartesian_deriv(t,1);

    %Calculate pseudo inverse of J to find angular rates
    J = Jacobian(q);
    JT = transpose(J);
    J_pseudo = TJ/(J*JT); %placeholder, this minimizes joint velocity
    
    %Calculate q_dot using pseudo inverse of J
    q_dot_1 = J_pseudo*xy_dot;
    
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