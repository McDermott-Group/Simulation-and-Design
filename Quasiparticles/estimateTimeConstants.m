function [tau_p, err_p, tau_r, err_r] = estimateTimeConstants(t, n_qp, plot_flag)

if ~exist('plot_flag', 'var')
    plot_flag = false;
end

% Extract poisoning time constant.
t_p = t(t < 0) - min(t);
n_p = n_qp(t < 0);
n_p = n_p - n_p(1);

n = @(t) interp1(t_p, n_p, t);

tau_p = fminsearch(@(t) abs(n(t) - max(n_p) * (1 - exp(-1))), 0);
err_p = 0;

% Extract recovery time constant.
t_r = t(t > 0);
t_r = t_r - min(t_r);
n_r = n_qp(t > 0);

n = @(t) interp1(t_r, n_r, t);

tau_r = fminsearch(@(t) abs(n(t) - max(n_r) * exp(-1)), 1);
err_r = 0;
end