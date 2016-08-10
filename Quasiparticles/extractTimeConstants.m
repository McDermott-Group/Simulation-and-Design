function [tau_p, err_p, tau_r, err_r] = extractTimeConstants(t, n_qp, plot_flag)

if ~exist('plot_flag', 'var')
    plot_flag = false;
end

% Extract poisoning time constant.
t_p = t(t < 0) - min(t);
n_p = n_qp(t < 0);
n_p = n_p - n_p(1);

% Cut off the intial part of the curve to improve the fit.
cut_off = 0 * max(n_p);
t_p(n_p < cut_off) = [];
n_p(n_p < cut_off) = [];
t_p = t_p - min(t_p);
n_p = n_p - cut_off;

pos = find(n_p == min(n_p(n_p > max(n_p) / exp(1))));
t_start = t_p(pos(1));

f_p = fit(t_p(:), n_p(:), 'a * (1 - exp(-b * x))',...
    'StartPoint', [max(n_p) , 1 / t_start], 'TolFun', 1e-35);

tau_p = 1 / f_p.b;
ci_p = confint(f_p);
err_p = max([abs(1 / ci_p(1, 2) - 1 / f_p.b),...
             abs(1 / ci_p(2, 2) - 1 / f_p.b)]);

if plot_flag
    figure
    hold on
    plot(t_p, n_p, 'LineWidth', 2)
    plot(t_p, f_p.a * (1 - exp(-t_p / tau_p)), 'LineWidth', 2)
    xlabel('Time (\tau_0)', 'FontSize', 14)
    ylabel('n_{\rm qp} / n_{\rm cp}', 'FontSize', 14)
    title('Poisoning Fit')
    legend('simulation', 'fit')
    grid on
    grid minor
    axis tight
end

% Extract recovery time constant.
t_r = t(t > 0);
t_r = t_r - min(t_r);
n_r = n_qp(t > 0);

% Cut off the intial part of the curve to improve the fit.
% cut_off = .8 * max(n_r);
% t_r(n_r > cut_off) = [];
% n_r(n_r > cut_off) = [];
% t_r = t_r - min(t_r);

cut_off = 0;
n_r(t_r < cut_off) = [];
t_r(t_r < cut_off) = [];
t_r = t_r - min(t_r);

f_r = fit(t_r(:), n_r(:), 'a * exp(-b * x)',...
    'StartPoint', [max(n_r), f_p.b], 'TolFun', 1e-35);
    
tau_r = 1 / f_r.b;
ci_r = confint(f_r);
err_r = max([abs(1 / ci_r(1, 2) - 1 / f_r.b),...
             abs(1 / ci_r(2, 2) - 1 / f_r.b)]);

if plot_flag
    figure
    hold on
    plot(t_r, n_r, 'LineWidth', 2)
    plot(t_r, f_r.a * exp(-t_r / tau_r), 'LineWidth', 2)
    xlabel('Time (\tau_0)', 'FontSize', 14)
    ylabel('n_{\rm qp} / n_{\rm cp}', 'FontSize', 14)
    title('Recovery Fit')
    legend('simulation', 'fit')
    grid on
    grid minor
    axis tight
end

end