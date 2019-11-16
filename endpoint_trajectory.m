clear
close all

load('FootPlacement.mat')
x = 0.8*x/100; z = 0.8*z/100; 
T = 0:0.01:timex(end);

Xf = fit(timex, x,'smoothingspline');
Zf = fit(timez, z,'smoothingspline');

figure()
plot(Xf(T),Zf(T),'k','LineWidth',2)
hold on

SSE_vec = zeros(2,7);
LEGEND = cell(1,8);
LEGEND{1} = 'Original Data';

for degree = 2:8
    poly_n = ['poly' num2str(degree)];
    x_smooth = fit(timex, x, poly_n);
    z_smooth = fit(timez, z, poly_n);
    
    plot(x_smooth(T),z_smooth(T))
    LEGEND{degree} = poly_n;
    
    SSE_x = 0;
    SSE_z = 0;
    for t = T
        SSE_x = SSE_x + (x_smooth(t) - Xf(t))^2;
        SSE_z = SSE_z + (z_smooth(t) - Zf(t))^2;
    end
    SSE_vec(1,degree - 1) = SSE_x;
    SSE_vec(2,degree - 1) = SSE_z;
end
legend(LEGEND)

figure(2)
plot(2:8,SSE_vec)
legend('SSE: X','SSE: Y')

