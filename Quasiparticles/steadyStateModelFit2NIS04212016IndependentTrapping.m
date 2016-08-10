function steadyStateModelFit2NIS04212016IndependentTrapping
%steadyStateModelFit2NIS04212016IndependentTrapping Fitting to
% the `NIS04212016` dataset using the two-region steady-state quasi-0D
% model.

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon_n and r_phonon_f dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 1.679e-04; r_phonon_n = 3.044e-01; r_phonon_f = 7.463e-02; c_qp = 7.463e-02; c_n = 1.705e-01; c_f = 3.499e-02; vol = 2.600e+04;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 250;

% Data.
data = load('NIS04212016.mat');
E_p_n = data.NearTrapPoisoning(:, 2);
nqp_p_n = data.NearTrapPoisoning(:, 4);

E_r_n = data.NearTrapRecovery(:, 2);
nqp_r_n = data.NearTrapRecovery(:, 4);

V_n = [E_p_n; E_r_n];
nqp_n = [nqp_p_n; nqp_r_n];

E_p_f = data.FarTrapPoisoning(:, 2);
nqp_p_f = data.FarTrapPoisoning(:, 4);

E_r_f = data.FarTrapRecovery(:, 2);
nqp_r_f = data.FarTrapRecovery(:, 4);

V_f = [E_p_f; E_r_f];
nqp_f = [nqp_p_f; nqp_r_f];

V_sim_n = V_n;
V_sim_f = V_f;
nqp_sim_n = NaN(size(V_sim_n));
nqp_sim_f = NaN(size(V_sim_f));
parfor k = 1:length(V_sim_n)
    if V_sim_n(k) > 1
        [~, ~, ~, ~, nqp] =...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V_sim_n(k), r_direct, r_phonon_n, c_qp, c_n, vol, N);
        nqp_sim_n(k) = max(nqp);
    else
        nqp_sim_n(k) = 0;
    end
    fprintf('*')
end
fprintf('\n')
parfor k = 1:length(V_sim_f)
    if V_sim_f(k) > 1
        [~, ~, ~, ~, nqp] = ...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V_sim_f(k), r_direct, r_phonon_f, c_qp, c_f, vol, N);
        nqp_sim_f(k) = max(nqp);
    else
        nqp_sim_f(k) = 0;
    end
    fprintf('*')
end
fprintf('\n')

h = figure;
hold on
plot(V_n, nqp_n, 'b^', V_f, nqp_f, 'bp',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'b')
plot(V_sim_n, nqp_sim_n, 'k^', V_sim_f, nqp_sim_f, 'kp',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'k')
hold off
xlabel('Normalized Injection Bias eV/\Delta', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'experiment: near resonator',...
        'experiment: far resonator',...
       ['simulation: r_{ph} = ', num2str(r_phonon_n, '%.3f'),...
                  '; c_{ph} = ', num2str(c_n, '%.3f')],...
       ['simulation: r_{ph} = ', num2str(r_phonon_f, '%.3f'),...
                  '; c_{ph} = ', num2str(c_f, '%.3f')]},...
        'Location', 'SouthEast')
title({'With Traps',...
    ['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0, ',...
     'c = ', num2str(c_qp, '%.3f'), ', ',...
     'vol = ', num2str(vol, '%.2e'), ' \mu{m}^3']})
axis tight
grid on
set(gca, 'yscale', 'Log')
set(gca, 'box', 'on')
savePDF(h, 'SimNearFarIndependentTrapping.pdf')

end