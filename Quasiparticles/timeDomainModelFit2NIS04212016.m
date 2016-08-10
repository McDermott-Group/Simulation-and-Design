function timeDomainModelFit2NIS04212016
%timeDomainModelFit2NIS04212016 Fitting to the FinalNIS04212016forSimulations
% dataset using the two-point time domain quasi-0D model.

data = load('NIS04212016.mat');

% V_p_f = data.FarTrapPoisoning(:, 1);
% E_p_f = data.FarTrapPoisoning(:, 2);
% tau_p_f = data.FarTrapPoisoning(:, 3);
% nqp_p_f = data.FarTrapPoisoning(:, 4);

% V_r_f = data.FarTrapRecovery(:, 1);
% E_r_f = data.FarTrapRecovery(:, 2);
% tau_r_f = data.FarTrapRecovery(:, 3);
% nqp_r_f = data.FarTrapRecovery(:, 4);

% V_p_n = data.NearTrapPoisoning(:, 1);
E_p_n = data.NearTrapPoisoning(:, 2);
tau_p_n = data.NearTrapPoisoning(:, 3);
nqp_p_n = data.NearTrapPoisoning(:, 4);

% V_r_n = data.NearTrapRecovery(:, 1);
E_r_n = data.NearTrapRecovery(:, 2);
tau_r_n = data.NearTrapRecovery(:, 3);
nqp_r_n = data.NearTrapRecovery(:, 4);

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 1.481e-05; r_phonon = 6.618e-01; c = 4.608e-02; vol = 2.600e+04;

V = [E_p_n; E_r_n; 3.5; 4.2; 6.7]; % in units of \Delta
Tph = 0.051; % K
tspan = [-500, 2000]; % in units of \tau_0

% Number of energy bins.
N = 150;

tau_p = NaN(size(V));
err_p = NaN(size(V));
tau_r = NaN(size(V));
err_r = NaN(size(V));
nqp = NaN(size(V));
parfor k = 1:length(V)
    [t, ~, ~, ~, n_qp] = ...
        twoRegionTimeDomainModelOptimized(Tph, tspan, V(k),...
        r_direct, r_phonon, c, vol, N);
    [tau_p(k), err_p(k), tau_r(k), err_r(k)] = ...
        estimateTimeConstants(t, n_qp, false);
    nqp(k) = max(n_qp);
    fprintf('*')
end
fprintf('\n')

tau0 = .438; % us, \tau_0 for aluminum from S. B. Kaplan et al.,
              % Phys. Rev. B 14, 4854 (1976)
tau_p = tau0 * tau_p;
tau_r = tau0 * tau_r;

F = median(tau_p_n ./ tau_p(1:length(tau_p_n)));

tau_p = F * tau_p;
tau_r = F * tau_r;

% load('SimNIS24062016.mat')

h = figure;
hold on
plot(E_p_n, tau_p_n, 'r^', 'MarkerSize', 9, 'LineWidth', 1.5)
plot(E_r_n, tau_r_n, 'r^', 'MarkerSize', 9, 'LineWidth', 1.5)
errorbar(V, tau_p, err_p, '*', 'MarkerSize', 10, 'LineWidth', 1.5)
errorbar(V, tau_r, err_r, '*', 'MarkerSize', 10, 'LineWidth', 1.5)
set(gca, 'xscale', 'Log')
hold off
xlabel('Normalized Injection Bias eV/\Delta', 'FontSize', 14)
ylabel('Time Constant (\mu{s})', 'FontSize', 14)
legend({'poisoning, near trap', 'recovery, near trap',...
    'poisoning simulation', 'recovery simulation'},...
    'Location', 'SouthWest')
title({'Time Constants',...
    ['r_{qp} = ', num2str(r_direct, '%.3e'), '/\tau_0',...
    ', r_{ph} = ', num2str(r_phonon, '%.3e'),...
    ', c = ', num2str(c, '%.3e'),...
    ', F = ', num2str(F, '%.3f')]})
axis tight
grid on
set(gca, 'box', 'on')
savePDF(h, 'SimNIS24062016_tau.pdf')
savefig(h, 'SimNIS240062016_tau.fig')

h = figure;
hold on
plot(E_p_n, tau_p_n, 'r^', 'MarkerSize', 9, 'LineWidth', 1.5)
plot(V, tau_p, 'mo', 'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'm')
set(gca, 'xscale', 'Log')
hold off
xlabel('Normalized Injection Bias eV/\Delta', 'FontSize', 14)
ylabel('Time Constant (\mu{s})', 'FontSize', 14)
legend({'poisoning, near trap', 'poisoning simulation'},...
    'Location', 'SouthWest')
title({'Time Constants',...
    ['r_{qp} = ', num2str(r_direct, '%.3e'), '/\tau_0',...
    ', r_{ph} = ', num2str(r_phonon, '%.3e'),...
    ', c = ', num2str(c, '%.3e'),...
    ', F = ', num2str(F, '%.3f')]})
axis tight
grid on
set(gca, 'box', 'on')
savePDF(h, 'SimNIS24062016_taup.pdf')
savefig(h, 'SimNIS240062016_taup.fig')
 
h = figure;
hold on
plot(E_p_n, nqp_p_n, 'r^', 'MarkerSize', 9, 'LineWidth', 1.5)
plot(E_r_n, nqp_r_n, 'r^', 'MarkerSize', 9, 'LineWidth', 1.5)
plot(V, nqp, 'mo', 'MarkerSize', 10, 'LineWidth', 1.5, 'MarkerFaceColor', 'm')
hold off
set(gca, 'xscale', 'Log', 'yscale', 'Log')
xlabel('Normalized Injection Bias eV/\Delta', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'poisoning, near trap', 'recovery, near trap', 'simulation'},...
    'Location', 'NorthWest')
title({'Quasiparticle Steady-State Density',...
    ['r_{qp} = ', num2str(r_direct, '%.3e'), '/\tau_0',...
    ', r_{ph} = ', num2str(r_phonon, '%.3e'),...
    ', c = ', num2str(c, '%.3e'),...
    ', F = ', num2str(F, '%.3f')]})
xlim([1 100])
axis tight
grid on
set(gca, 'box', 'on')
savePDF(h, 'SimNIS24062016_nqp.pdf')
savefig(h, 'SimNIS24062016_nqp.fig')

save('SimNIS24062016.mat', 'E_p_n', 'tau_p_n', 'nqp_p_n',...
    'E_r_n', 'tau_r_n', 'nqp_r_n', 'V', 'nqp', 'tau_p', 'tau_r', 'F',...
    'err_p', 'err_r')

end