function plotSteadyStateDistributionsDifferentBias
%plotSteadyStateDistributionsDifferentBias Quasiparticle distribution plots
% for different biases.

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon dimensionless
% c dimensionless
% vol in units of um^3
r_direct =1e-05; r_phonon = 5e-01; c = 0; vol = 2.6e+04; % um^3

N = 250;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

V = [1.25, 2.5, 5, 10];

e = cell(length(V));
f_res = cell(length(V));
f_nis = cell(length(V));
parfor k = 1:length(V)
    [~, e{k}, ~, f_res{k}, ~, ~, ~, ~, ~, ~, ~, ~, f_nis{k}] = ...
        twoRegionSteadyStateModelOptimized(Tph, tspan, V(k),...
        r_direct, r_phonon, c, c, vol, N);
end

h = figure;
hold on
for k = 1:length(V)
    fnis = f_nis{k};
    plot(e{k}, fnis(end, :), 'LineWidth', 3)
end
ax = gca;
ax.ColorOrderIndex = 1;
for k = 1:length(V)
    fres = f_res{k};
    plot(e{k}, fres(end, :), '--', 'LineWidth', 3)
end
hold off
xlabel('Energy (\epsilon=E/\Delta)', 'FontSize', 14)
ylabel('Occupation Numbers f', 'FontSize', 14)
legends = cell(0);
for k = 1:length(V)
    legends{k} = ['NIS: V = ', num2str(V(k)), ' \Delta/e'];
end
for k = 1:length(V)
    legends{length(V) + k} = ['resonator: V = ', num2str(V(k)), ' \Delta/e'];
end
legend(legends)
set(gca, 'yscale', 'Log')
axis tight
ylim([1e-6, 2e-3])
xlim([1, max(V) / 4])
grid on

savePDF(h, ['SimDistributions_c', strrep(num2str(c), '.', 'p'), '.pdf'])
end