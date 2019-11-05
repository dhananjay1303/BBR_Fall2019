%% Attempt at bootleg numerical differentiation

i = 50;
t = 0:0.01:1;
y = exp(t);

POI = y(i:i+6);
TOI = t(i:i+6);

fitObj = fit(TOI',POI','poly5');
[fitObjPrime, fitObjDoublePrime] = differentiate(fitObj, TOI');

disp(num2str(fitObjPrime(end)));
disp(num2str(cos(t(i+4))));

disp(' ');

disp(num2str(fitObjDoublePrime(end)));
disp(num2str(-sin(t(i+6))));