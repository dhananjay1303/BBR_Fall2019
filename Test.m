clear
close all

load('FootPlacement.mat')
Xf = fit(timex, x/100,'smoothingspline'); %poly8
Zf = fit(timez, z/100,'smoothingspline'); %poly4

X = [0];
Z = [0];
t = [0];

while t(end) < timex(end)
    X = [X Xf(t(end))];
    Z = [Z Zf(t(end))];
    t = [t t(end)+0.01];
end
X(1) = [];
Z(1) = [];

plot(X,Z)

figure()
plot(t, differentiate(Xf, t), t, differentiate(Zf, t))