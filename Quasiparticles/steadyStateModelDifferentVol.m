function steadyStateModelDifferentVol
%steadyStateModelDifferentVol Explore the quasiparticle steady-state
% density at different injection volumes.

r_direct = 1e-05; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon = 1; % dimensionless
c = 0; % dimensionless
vol = [0.1 1 10 100] * 1.000e+03; % um^3
V = linspace(1.001, 30, 30); % in units of \delt

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

nqp = NaN(length(V), length(vol));
P = NaN(size(nqp));
for krqp = 1:length(vol)
    volsim = vol(krqp);
    parfor kV = 1:length(V)
        Vsim = V(kV);
        [~, ~, ~, ~, n_qp, ~, P(kV, krqp)] =...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            Vsim, r_direct, r_phonon, c, c, volsim, N);
        nqp(kV, krqp) = max(n_qp);
        fprintf('*')
    end
end

h = figure;
hold on
for k = 1:length(vol)
    plot(P(:, k), nqp(:, k), 'MarkerSize', 10, 'LineWidth', 2)
end
xlabel('Power (W)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legends = cell(0);
for k = 1:length(vol)
    legends{k} = ['v = ', num2str(vol(k), '%.2e'), ' \mu{m}^3'];
end
legend(legends, 'Location', 'SouthEast')
title(['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0', ', ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f')])
axis tight
set(gca, 'xscale', 'Log')
set(gca, 'yscale', 'Log')
grid on
set(gca, 'box', 'on')
savePDF(h, 'v.pdf')
end