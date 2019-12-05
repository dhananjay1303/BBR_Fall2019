clear
close all

%% Import Experimental Data
load('FootPlacement.mat')
Xf = fit(timex, 0.8*(x/100),'smoothingspline');
Zf = fit(timez, 0.8*(z/100),'smoothingspline');

%% Define MPC Object and Parameters
t = linspace(0,timex(end),62); 

%-- Parameter Definitions --%
Ts = t(2) - t(1);           % Sample Time
last_state_dot = 0;         % Record last state derivative

%Leg Geometry
lh = 0.424; lk = 0.425; la = 0.034;
leg_length = struct('hip',lh,'knee',lk,'ankle',la);
q_star = deg2rad([0;10;90]); %comfortable positions

%Parameter Cell
params = nlmpcmoveopt; %Non-linear mpc "options" object
params.Parameters = {Ts, last_state_dot, leg_length, q_star};
n_params = length(params.Parameters);

%MPC Object
LEG = nlmpc(3,2,3); %3 states, 2 outputs, 3 inputs

%Model Definition
LEG.Ts = Ts;
LEG.ControlHorizon = 8;
LEG.PredictionHorizon = 8;
LEG.Model = struct('StateFcn',@myStateFunction,'OutputFcn',@myOutputFunction,'NumberOfParameters',n_params);
LEG.Jacobian = struct('OutputFcn',@Jacobian);
LEG.Optimization = struct('CustomCostFcn',@myCostFunction,'ReplaceStandardCost',false);
LEG.Weights.OutputVariables = [40 40];
LEG.Weights.ManipulatedVariablesRate = [0.01 0.01 0.01];

%'CustomCostFcn',@myCostFunction,

%Joint Limits
LEG.States(1).Min = deg2rad(-45);
LEG.States(1).Max = deg2rad( 45);
LEG.States(2).Min = deg2rad(  0);
LEG.States(2).Max = deg2rad( 80);
LEG.States(3).Min = deg2rad( 60);
LEG.States(3).Max = deg2rad(170);

%% Simulation 
N = length(t);

%Initial State
state = deg2rad([-5;10;75]);
qh_0 = state(1);
qk_0 = state(2);
qa_0 = state(3);

%Initial End Point Position
x0 = lh*sin(qh_0) + lk*sin(qh_0-qk_0) + la*sin(qh_0-qk_0+qa_0);
z0 = -lh*cos(qh_0) - lk*cos(qh_0-qk_0) - la*cos(qh_0-qk_0+qa_0);

%Moving displacement vector to intertial frame
delta_x = x0 - Xf(0);
delta_z = z0 - Zf(0);

%Live plot initialization
figure()
plot(Xf(t)+delta_x,Zf(t)+delta_z)
hold on
axis equal

lastmv = zeros(3,1);                    % Controller input from last time step
state_vector = zeros(3,N);              % Vector for recording state time history
step_count = LEG.PredictionHorizon - 1; % How far to look into the future

%-- Loop through euler integration simulation --%
for i = 1:N
    
    %Set reference trajectory until the prediction horizon
    try
        ref_point = [Xf(t(i:i+step_count))+delta_x Zf(t(i:i+step_count))+delta_z];
    catch
        ref_point = [Xf(t(i:end))+delta_x Zf(t(i:end))+delta_z];
    end
    
    %Decide optimal control input
    mv = nlmpcmove(LEG,state,lastmv,ref_point,[],params);
    
    %Euler Integration using Control Input
    state_dot = myStateFunction(0,mv,Ts,last_state_dot);
    state = state + Ts*state_dot;
    lastmv = mv;
    last_state_dot = state_dot;
    
    %Live plot
    state_vector(:,i) = state';
    output = myOutputFunction(state,0,0,0,leg_length,0);
    plot(ref_point(1,1),ref_point(1,2),'cx')
    plot(output(1),output(2),'gd')
    drawnow
end

%% Plot Results

for i = 1:length(state_vector)
    logic_array = [(state_vector(:,i) ~= zeros(3,1)); mod(i,2)==0];
    if logic_array
        plotLeg(state_vector(:,i))
    end
end

save('SimulationResults.mat','state_vector','-append')
%% Functions

function state_dot = myStateFunction(~,u,~,~,~,~)
    
    %We define the input u to be the a column vector containing the angular
    %velocity of all 3 joint angles
    state_dot = u;
    
end

function output = myOutputFunction(x,~,~,~,leg_length,~)

    qh = x(1); qk = x(2); qa = x(3);
    
    %Leg geometry
    lh = leg_length.hip;
    lk = leg_length.knee;
    la = leg_length.ankle;

    x0 = lh*sin(qh) + lk*sin(qh-qk) + la*sin(qh-qk+qa);
    z0 = -lh*cos(qh) - lk*cos(qh-qk) - la*cos(qh-qk+qa);

    output = [x0; z0];

end

function J = Jacobian(x,~,~,~,leg_length,~)

    qh = x(1); qk = x(2); qa = x(3);
    
    %Leg geometry
    lh = leg_length.hip;
    lk = leg_length.knee;
    la = leg_length.ankle;
    
    %Jacobian Definition
    J = zeros(2,3);
    J(1,1) = lh*cos(qh) + lk*cos(qh-qk) + la*cos(qh-qk+qa);
    J(1,2) = -lk*cos(qh-qk) - la*cos(qh-qk+qa);
    J(1,3) = la*cos(qh-qk+qa);
    J(2,1) = lh*sin(qh) + lh*sin(qh-qk) + la*sin(qh-qk+qa);
    J(2,2) = -lh*sin(qh-qk) - la*sin(qh-qk+qa);
    J(2,3) = la*sin(qh-qk+qa);
end

function cost = myCostFunction(X,U,~,data,Ts,~,~,q_star)
    p = data.PredictionHorizon;

    %Minimize Jerk
    vel_h = U(2:p+1,1);
    vel_k = U(2:p+1,2);
    vel_a = U(2:p+1,3);
    
    acc_h = diff(vel_h)/Ts;
    acc_k = diff(vel_k)/Ts;
    acc_a = diff(vel_a)/Ts;
    
    jerk_h = diff(acc_h)/Ts;
    jerk_k = diff(acc_k)/Ts;
    jerk_a = diff(acc_a)/Ts;
    
    jerk_cost = sum(jerk_h.^2) + sum(jerk_k.^2) + sum(jerk_a.^2);
    
    %Scaling to be on the same order as the other cost functions
    %(This was somewhat arbitrary) 
    jerk_cost = jerk_cost/(1000*(p+1)); 
    
    %Avoid uncomfortable positions
    qh = X(2:p+1,1); qk = X(2:p+1,2); qa = X(2:p+1,3);
    
    hip_cost = sum((qh - q_star(1)).^2);
    knee_cost = sum((qk - q_star(2)).^2);
    ankle_cost = 5*sum((qa - q_star(3)).^2);
    
    %Find Total Cost
    cost = jerk_cost + hip_cost + knee_cost + ankle_cost;
%     cost = cost;

end

function plotLeg(q)

    lh = 0.424; lk = 0.425; la = 0.034;
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
    %drawnow
end