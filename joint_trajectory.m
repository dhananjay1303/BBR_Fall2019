% close all
clear x z T Xf Zf

load('FootPlacement.mat')
x = 0.8*x/100; z = 0.8*z/100; 
T = 0:0.001:timex(end);

Xf = fit(timex, x,'smoothingspline');
Zf = fit(timez, z,'smoothingspline');

figure()
plot(Xf(T),Zf(T)-0.92)
hold on

qh = [-0.087;0.08;0.145;0.3146;0.3259;0.3308];
qk = [0.1745;0.79;0.895;0.4295;0.1456;0.0685];
qa = [1.309;1.64;1.67;1.56;1.311;0.9407];
t = [0;0.1924;0.2388;0.468;0.5643;0.6188];

qh_fit = fit(t,qh,'poly5');
qk_fit = fit(t,qk,'poly5');
qa_fit = fit(t,qa,'poly5');

count = 0;
for time = T
    count = count + 1;
    if count == 3
        plotLeg([qh_fit(time);qk_fit(time);qa_fit(time)]);
        pause(0.05)
        count = 0;
    end
end

qh_dot_jerk = differentiate(qh_fit,T);
qk_dot_jerk = differentiate(qk_fit,T);
qa_dot_jerk = differentiate(qa_fit,T);

figure(3)
plot(T,(qh_dot_jerk.^2 + qk_dot_jerk.^2 + qa_dot_jerk.^2))
legend('Minimum Velocity','Minimum Kinetic Energy','Minimum Jerk')

jerk_cost_minJerk = plotJerk(qh_fit, qk_fit, qa_fit,T,1)/length(T);

%-- Plotting Leg Position --%
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