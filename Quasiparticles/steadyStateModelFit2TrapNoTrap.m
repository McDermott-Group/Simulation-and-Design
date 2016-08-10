function steadyStateModelFit2TrapNoTrap
%steadyStateModelFit2TrapNoTrap Fitting to the `TrapNoTrap` dataset using
% the two-region steady-state quasi-0D model.

r_direct = 1.501e-05; r_phonon = 1.076e+00; c = 1.159e-02; vol = 2.600e+04;
r_direct_no_tr = r_direct; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon_no_tr = r_phonon; % dimensionless
c_no_tr = c; % dimensionless
vol_no_tr = vol; % um^3

r_direct = 3.963e-05; r_phonon = 8.102e-01; c = 8.885e-02; vol = 2.600e+04;
r_direct_tr = r_direct; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon_tr = r_phonon; % dimensionless
c_tr = c; % dimensionless
vol_tr = vol; % um^3

delta = 0.18e-3; % eV (aluminum superconducting gap)
Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 250;

data = load('TrapNoTrap.mat');

V_tr = data.Trap(:, 5) / delta;
P_tr = data.Trap(:, 6);
nqp_tr = data.Trap(:, 8) - min(data.Trap(:, 8));

V_no_tr = data.NoTrap(:, 5) / delta;
P_no_tr = data.NoTrap(:, 6);
nqp_no_tr = data.NoTrap(:, 8) - min(data.NoTrap(:, 8));

figure
hold on
plot(V_no_tr, nqp_no_tr, 'k^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'k')
plot(V_tr, nqp_tr, 'r^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'r')
hold off
set(gca, 'xscale', 'Log', 'yscale', 'Log')
xlabel('Normalized Injection Bias eV/\Delta', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend('no traps', 'with traps', 'Location', 'SouthEast')
title('Experiment')
axis tight
grid on
set(gca, 'box', 'on')

figure
hold on
plot(P_no_tr, nqp_no_tr, 'k^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'k')
plot(P_tr, nqp_tr, 'r^',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'r')
hold off
set(gca, 'xscale', 'Log', 'yscale', 'Log')
xlabel('Power (W)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend('no traps', 'with traps', 'Location', 'SouthEast')
title('Experiment')
axis tight
grid on
set(gca, 'box', 'on')

V_sim_no_tr = [1.005; V_no_tr; 30.5];
V_sim_tr = [1.005; V_tr; 30.5];
nqp_sim_no_tr = NaN(size(V_sim_no_tr));
nqp_sim_tr = NaN(size(V_sim_tr));
P_sim_no_tr = NaN(size(V_sim_no_tr));
P_sim_tr = NaN(size(V_sim_tr));
Gamma_tr_nis_no_tr = NaN(size(V_sim_no_tr));
Gamma_tr_res_no_tr = NaN(size(V_sim_no_tr));
Gamma_tr_nis_tr = NaN(size(V_sim_tr));
Gamma_tr_res_tr = NaN(size(V_sim_tr));
Gamma_r_nis_no_tr = NaN(size(V_sim_no_tr));
Gamma_r_res_no_tr = NaN(size(V_sim_no_tr));
Gamma_r_nis_tr = NaN(size(V_sim_tr));
Gamma_r_res_tr = NaN(size(V_sim_tr));
parfor k = 1:length(V_sim_no_tr)
    if V_sim_no_tr(k) > 1
        [~, ~, ~, ~, nqp, ~, P_sim_no_tr(k),...
            Gamma_tr_nis_no_tr(k), Gamma_tr_res_no_tr(k),...
            Gamma_r_nis_no_tr(k), Gamma_r_res_no_tr(k)] =...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V_sim_no_tr(k), r_direct_no_tr, r_phonon_no_tr, c_no_tr,...
            c_no_tr, vol_no_tr, N);
        nqp_sim_no_tr(k) = max(nqp);
    else
        nqp_sim_no_tr(k) = 0;
    end
    fprintf('*')
end
fprintf('\n')
parfor k = 1:length(V_sim_tr)
    if V_sim_tr(k) > 1
        [~, ~, ~, ~, nqp, ~, P_sim_tr(k),...
            Gamma_tr_nis_tr(k), Gamma_tr_res_tr(k),...
            Gamma_r_nis_tr(k), Gamma_r_res_tr(k)] = ...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V_sim_tr(k), r_direct_tr, r_phonon_tr, c_tr,...
            c_tr, vol_tr, N);
        nqp_sim_tr(k) = max(nqp);
    else
        nqp_sim_tr(k) = 0;
    end
    fprintf('*')
end
fprintf('\n')

figure
semilogy(V_no_tr, nqp_no_tr, 'k^', V_sim_no_tr, nqp_sim_no_tr, 'm*',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'k')
xlabel('Normalized Injection Bias eV/\Delta', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'no traps, experiment',...
    ['r_{qp} = ', num2str(r_direct_no_tr, '%.2e'), '/\tau_0, ',...
    'r_{ph} = ', num2str(r_phonon_no_tr, '%.3f'), ', ',...
    'c = ', num2str(c_no_tr, '%.3f'), ', ',...
    'vol = ', num2str(vol_no_tr, '%.2e'), ' \mu{m}^3']},...
    'Location', 'SouthEast')
title('No Traps')
axis tight
grid on
 
figure
semilogy(V_tr, nqp_tr, 'r^', V_sim_tr, nqp_sim_tr, 'm*',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'r')
xlabel('Normalized Injection Bias eV/\Delta', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'no trap, experiment',...
    ['r_{qp} = ', num2str(r_direct_tr, '%.2e'), '/\tau_0, ',...
    'r_{ph} = ', num2str(r_phonon_tr, '%.3f'), ', ',...
    'c = ', num2str(c_tr, '%.3f'), ', ',...
    'vol = ', num2str(vol_tr, '%.2e'), ' \mu{m}^3']},...
    'Location', 'SouthEast')
title('With Traps')
axis tight
grid on

h = figure;
loglog(P_no_tr, nqp_no_tr, 'k^', P_sim_no_tr, nqp_sim_no_tr, 'm*',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'k')
xlabel('Power (W)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'no trap, experiment',...
    ['r_{qp} = ', num2str(r_direct_no_tr, '%.1e'), '/\tau_0, ',...
    'r_{ph} = ', num2str(r_phonon_no_tr, '%.2f'), ', ',...
    'c = ', num2str(c_no_tr, '%.2f'), ', ',...
    'vol = ', num2str(vol_no_tr, '%.1e'), ' \mu{m}^3']},...
    'Location', 'SouthEast')
title('No Traps')
axis tight
grid on
set(gca, 'box', 'on')
savePDF(h, 'SimNoTrap.pdf')
save('SimNoTrap.mat', 'V_no_tr', 'P_no_tr', 'nqp_no_tr',...
    'V_sim_no_tr', 'P_sim_no_tr', 'nqp_sim_no_tr')

h = figure;
loglog(P_tr, nqp_tr, 'r^', P_sim_tr, nqp_sim_tr, 'm*',...
    'MarkerSize', 9, 'LineWidth', 2, 'MarkerFaceColor', 'r')
xlabel('Power (W)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legend({'trap, experiment',...
    ['r_{qp} = ', num2str(r_direct_tr, '%.1e'), '/\tau_0, ',...
    'r_{ph} = ', num2str(r_phonon_tr, '%.2f'), ', ',...
    'c = ', num2str(c_tr, '%.2f'), ', ',...
    'vol = ', num2str(vol_tr, '%.1e'), ' \mu{m}^3']},...
    'Location', 'SouthEast')
title('With Traps')
axis tight
grid on
set(gca, 'box', 'on')
savePDF(h, 'SimTrap.pdf')
save('SimTrap.mat', 'V_tr', 'P_tr', 'nqp_tr',...
    'V_sim_tr', 'P_sim_tr', 'nqp_sim_tr')

h = figure;
loglog(P_sim_no_tr, Gamma_tr_nis_no_tr, 'k^',...
       P_sim_no_tr, Gamma_tr_res_no_tr, 'b^',...
       P_sim_tr, Gamma_tr_nis_tr, 'r^',...
       P_sim_tr, Gamma_tr_res_tr, 'm^',...
    'MarkerSize', 9, 'LineWidth', 2)
xlabel('Power (W)', 'FontSize', 14)
ylabel('Trapping Rate \Gamma_{tr} (s^{-1})', 'FontSize', 14)
legend({'no traps, NIS', 'no traps, resonator',...
        'with traps, NIS', 'with traps, resonator'},...
        'Location', 'NorthWest')
title('Trapping Rate')
axis tight
grid on
set(gca, 'box', 'on')
savePDF(h, 'SimGammaTrapping.pdf')

h = figure;
loglog(P_sim_no_tr, Gamma_r_nis_no_tr, 'k^',...
       P_sim_no_tr, Gamma_r_res_no_tr, 'b^',...
       P_sim_tr, Gamma_r_nis_tr, 'r^',...
       P_sim_tr, Gamma_r_res_tr, 'm^',...
    'MarkerSize', 9, 'LineWidth', 2)
xlabel('Power (W)', 'FontSize', 14)
ylabel('Recombination Rate \Gamma_{rec} (s^{-1})', 'FontSize', 14)
legend({'no traps, NIS', 'no traps, resonator',...
        'with traps, NIS', 'with traps, resonator'},...
        'Location', 'SouthEast')
title('Recombination Rate')
axis tight
grid on
set(gca, 'box', 'on')
savePDF(h, 'SimGammaRecombination.pdf')

end