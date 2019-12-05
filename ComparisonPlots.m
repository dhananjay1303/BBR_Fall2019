clear
close all

load('JointData.mat')
load('SimulationResults.mat')
load('FootPlacement.mat')

x = linspace(52,100,62);
t = linspace(0,timex(end),62);

%% Hip Plot

figure()
hold on
grid on
plot(t,deg2rad(HipFit(x)*0.47/0.71),'-')
plot(TOUT_JOINT_VEL,JOINTSPACE_VEL(:,1),'x-')
plot(TOUT_JOINT_KE,JOINTSPACE_KE(:,1),'d-')
plot(t,qh_fit(t),'o-')
plot(t,state_vector(1,:),'s-')
legend('Experimental','Min Velocity','Min Kinetic Energy','Min Jerk (Analytic)','Min Jerk (MPC)','Location','SouthEast')
title('Hip Angle (Scaled by Stride Length)')
ylabel('Radians')
xlabel('Seconds')

%% Knee Plot
knee_offset = 0;

figure()
hold on
grid on
plot(t,deg2rad(KneeFit(x)/1.25-knee_offset),'-')
plot(TOUT_JOINT_VEL,JOINTSPACE_VEL(:,2),'x-')
plot(TOUT_JOINT_KE,JOINTSPACE_KE(:,2),'d-')
plot(t,qk_fit(t),'o-')
plot(t,state_vector(2,:),'s-')
legend('Experimental','Min Velocity','Min Kinetic Energy','Min Jerk (Analytic)','Min Jerk (MPC)')
title('Knee Angle (Scaled by Stride Length)')
ylabel('Radians')
xlabel('Seconds')

%% Ankel Plot
ankle_offset = -90;

figure()
hold on
grid on
plot(t,deg2rad(-AnkleFit(x)-ankle_offset),'-')
plot(TOUT_JOINT_VEL,JOINTSPACE_VEL(:,3),'x-')
plot(TOUT_JOINT_KE,JOINTSPACE_KE(:,3),'d-')
plot(t,qa_fit(t),'o-')
plot(t,state_vector(3,:),'s-')
legend('Experimental','Min Velocity','Min Kinetic Energy','Min Jerk (Analytic)','Min Jerk (MPC)')
title('Ankle Angle')
ylabel('Radians')
xlabel('Seconds')

%% Data Plot
figure()
hold on
for i = 1:length(x)
    plotLeg(deg2rad([HipFit(x(i)) KneeFit(x(i)) -AnkleFit(x(i))+90]),(i/80)*[1 1 1])
end

%% Plotting Leg Position

function plotLeg(q,color)

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

    plot([xh xk xa x],[yh yk ya y],'color',color)
    axis equal
    %drawnow
end