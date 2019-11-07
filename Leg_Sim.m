clear
close all

%% Import Experimental Data
global Xf Zf
load('FootPlacement.mat')
Xf = fit(timex, x/100,'smoothingspline');
Zf = fit(timez, z/100,'smoothingspline');

%% System Parameter Definitions
global lh lk la qstar count bigCount

%Lengths of Leg Sections (not set yet)
lh = 0.8; %Upper leg length               (m)
lk = 0.6; %Lower leg length               (m)
la = 0.1; %Foot length from ankle to cg   (m)

%Joint Limits for Comfort
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

figure()
hold on

[TOUT_END, ENDPOINT] = ode45(@v, [0 timex(end)], [x0;z0]);
plot(ENDPOINT(:,1),ENDPOINT(:,2))
%plot(Xf(timex),Zf(timex))
drawnow

[TOUT_JOINT, JOINTSPACE] = ode45(@state_deriv,[0 timex(end)],[qh_0;qk_0;qa_0]);
%ode45(@state_deriv,[0 timex(end)],[qh_0;qk_0;qa_0]);

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

function velocity = v(t,~)
    global Xf Zf
    
    velocity = [differentiate(Xf,t); differentiate(Zf,t)];
end

function q_dot = state_deriv(t,q)
    
    global qstar Xf Zf count bigCount

    xz_dot = [differentiate(Xf,t);differentiate(Zf,t)];
    
    %Calculate pseudo inverse of J to find angular rates
    J = Jacobian(q);
    JT = transpose(J);
    J_pseudo = JT/(J*JT);
    
    q_dot_1 = J_pseudo*xz_dot;
    
    %-- Optimize to avoid uncomfortable joint positions --%
    X = 1; %Rate of convergence
    
    %Project solution onto null space
    n_null = null(J);
    q_dot_null = -X*dot((q-qstar)/norm(q-qstar), n_null)*n_null;
    
    %Determine the final q_dot solution
    q_dot = q_dot_1 + q_dot_null;
    
    %Only plot every 4th position
    if count ==4
        plotLeg(q)
        count = 0;
    else
        count = count +1;
    end
    
    bigCount = bigCount + 1;
    if bigCount > 500
        %assert(false)
    end
end

%% Plotting Arm Position

function plotLeg(q)

    global lh lk la bigCount
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

    plot([xh xk xa x],[yh yk ya y],'color',[1 1 1].*max(0.7-bigCount/100,0))
end