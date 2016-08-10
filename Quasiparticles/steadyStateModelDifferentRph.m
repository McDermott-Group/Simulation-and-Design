function steadyStateModelDifferentRph
%steadyStateModelDifferentRph Explore the quasiparticle steady-state
% density at different phonon fractions.

r_direct = 1e-05; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon = [.001, .01, 0.1, 1]; % dimensionless
c = 0; % dimensionless
vol = 2.6e+04; % um^3
V = linspace(1.001, 30, 30); % in units of \delta

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

nqp = NaN(length(V), length(r_phonon));
P = NaN(size(nqp));
for krqp = 1:length(r_phonon)
    rphsim = r_phonon(krqp);
    parfor kV = 1:length(V)
        Vsim = V(kV);
        [~, ~, ~, ~, n_qp, ~, P(kV, krqp)] =...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            Vsim, r_direct, rphsim, c, c, vol, N);
        nqp(kV, krqp) = max(n_qp);
        fprintf('*')
    end
end

h = figure;
hold on
for k = 1:length(r_phonon)
    plot(P(:, k), nqp(:, k), 'MarkerSize', 10, 'LineWidth', 2)
end
xlabel('Power (W)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legends = cell(0);
for k = 1:length(r_phonon)
    legends{k} = ['r_{ph} = ', num2str(r_phonon(k), '%.2e')];
end
legend(legends, 'Location', 'SouthEast')
title(['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0', ', ',...
     'c = ', num2str(c, '%.3f'), ', ',...
     'v = ', num2str(vol, '%.2e'), ' \mu{m}^3'])
axis tight
set(gca, 'xscale', 'Log')
set(gca, 'yscale', 'Log')
grid on
set(gca, 'box', 'on')
savePDF(h, 'rph.pdf')
end