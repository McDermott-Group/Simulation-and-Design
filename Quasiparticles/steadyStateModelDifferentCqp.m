function steadyStateModelDifferentCqp
%steadyStateModelDifferentCqp Explore the quasiparticle steady-state
% density at different trapping efficiency.

r_direct = 1e-05; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon = 1; % dimensionless
c = [0, .01, .05, .1]; % dimensionless
vol = 2.6e+04; % um^3
V = linspace(1.001, 30, 30); % in units of \delta

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

nqp = NaN(length(V), length(c));
P = NaN(size(nqp));
for krqp = 1:length(c)
    csim = c(krqp);
    parfor kV = 1:length(V)
        Vsim = V(kV);
        [~, ~, ~, ~, n_qp, ~, P(kV, krqp)] =...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            Vsim, r_direct, r_phonon, csim, 0, vol, N);
        nqp(kV, krqp) = max(n_qp);
        fprintf('*')
    end
end

h = figure;
hold on
for k = 1:length(c)
    plot(P(:, k), nqp(:, k), 'MarkerSize', 10, 'LineWidth', 2)
end
xlabel('Power (W)', 'FontSize', 14)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 14)
legends = cell(0);
for k = 1:length(c)
    legends{k} = ['c^{NIS} = ', num2str(c(k), '%.2e')];
end
legend(legends, 'Location', 'SouthEast')
title(['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0', ', ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c^{res} = ', num2str(0, '%.3f'), ', ',...
     'v = ', num2str(vol, '%.2e'), ' \mu{m}^3'])
axis tight
set(gca, 'xscale', 'Log')
set(gca, 'yscale', 'Log')
grid on
set(gca, 'box', 'on')
savePDF(h, 'cqp.pdf')
end