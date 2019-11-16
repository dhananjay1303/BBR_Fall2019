function jerk_cost = plotJerk(qh_fit,qk_fit,qa_fit,T,PLT)

%ensure T is a column vector
dims = size(T);
if dims(2) > 1
    T = transpose(T);
end

%Calculate Jerk Profile
[~, qh_acc] = differentiate(qh_fit,T);
[~, qk_acc] = differentiate(qk_fit,T);
[~, qa_acc] = differentiate(qa_fit,T);

%Fit acceleration to smooth curve fit
qh_acc_fit = fit(T,qh_acc,'smoothingspline');
qk_acc_fit = fit(T,qk_acc,'smoothingspline');
qa_acc_fit = fit(T,qa_acc,'smoothingspline');

%Find the derivative of acceleration
qh_jerk = differentiate(qh_acc_fit,T);
qk_jerk = differentiate(qk_acc_fit,T);
qa_jerk = differentiate(qa_acc_fit,T);

%Sum jerk cost over all 3 joints
jerk_cost = sum(qh_jerk.^2 + qk_jerk.^2 + qa_jerk.^2);

%Plot jerk cost over time
if PLT
    figure(6)
    semilogy(T,(qh_jerk.^2 + qk_jerk.^2 + qa_jerk.^2))
    hold on
    xlim([T(1) T(end)])
end
