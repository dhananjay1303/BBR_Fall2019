function poly_X = generateMinimumJerk()

% x,z = a0 + a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5

var = load('FootPlacement.mat');
X = var.x/100;
Z = var.z/100;
tx = var.timex;
tz = var.timez;
% Need 6 known points each to solve for the constants of the polynomial

%Initial Point
x0 = X(1);
z0 = Z(1);
t0 = 0;

%Furthest point the foot moves backward
x_min = min(X);
t_min_x = tx(X == x_min);
z_x_min = interp1(tz,Z,t_min_x);

%Highest Point (z)
z_max = max(Z);
t_max_z = tz(Z == z_max);
x_z_max = interp1(tx,X,t_max_z);

%Minimum Z after apex
apex = find(Z == z_max);
z_min_local = min(Z(apex:end));
t_min_local_z = tz(Z == z_min_local);
x_min_local = interp1(tx,X,t_min_local_z);

%Some other point
t_sample = 0.078;
x_sample = interp1(tx,X,t_sample);
z_sample = interp1(tz,Z,t_sample);

%Final Point
t_final_x = tx(end);
t_final_z = tz(end);
x_final = X(end);
z_final = Z(end);

%Find minimum jerk x-trajectory
X_MAT = [tVec(0);tVec(t_min_x);tVec(t_max_z);tVec(t_min_local_z);tVec(t_sample);tVec(t_final_x)];
RHS_X = [x0;x_min;x_z_max;x_min_local;x_sample;x_final];

poly_X = X_MAT\RHS_X;
x_vec = zeros(1,length(tx));

for i = 1:length(tx)
    x_vec(i) = polyval(poly_X,tx(i));
end
close all
figure(1)
plot(tx,x_vec,tx,X)

%Find minimum jerk z-trajectory
Z_MAT = [tVec(0);tVecPrime(t_max_z);tVec(t_max_z);tVec(t_min_local_z);tVecPrime(t_min_local_z);tVec(t_final_z)];
RHS_Z = [z0;0;z_max;z_min_local;0;z_final];

poly_Z = Z_MAT\RHS_Z;

tz_vec = 0:0.001:tx(end);
z_vec = zeros(1,length(tz_vec));

for i = 1:length(tz_vec)
    z_vec(i) = polyval(poly_Z,tz_vec(i));
end
figure(2)
plot(tz_vec,z_vec,tz,Z)

figure(3)

plot(polyval(poly_X,tx),polyval(poly_Z,tx))

    function timeVec = tVec(t)
        timeVec = [t^5 t^4 t^3 t^2 t 1];
    end

    function timeVecPrime = tVecPrime(t)
        timeVecPrime = [5*t^4 4*t^3 3*t^2 2*t 1 0];
    end
end
