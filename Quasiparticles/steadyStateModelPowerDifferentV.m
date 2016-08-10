function steadyStateModelPowerDifferentV
%steadyStateModelPowerDifferentV Explore the quasiparticle power budget
%at different injection voltages.

r_direct = 2.5e-05; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon = 0.9; % dimensionless
c = 0; % dimensionless
vol = 3.8e+04; % um^3
V = linspace(1.01, 10, 76); % in units of \delta

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

Prec = nan(size(V));
Psct = nan(size(V));
Psct2D = nan(size(V));
Ptrp = nan(size(V));
Ptrp2D = nan(size(V));
parfor kV = 1:length(V)
    Vsim = V(kV);
    [~, ~, ~, ~, ~, ~, ~, ...
        Prec(kV), Psct(kV), Psct2D(kV), Ptrp(kV), Ptrp2D(kV)] =...
        twoRegionSteadyStateModel(Tph, tspan,...
        Vsim, r_direct, r_phonon, c, vol, N, false);
    fprintf('*')
end

h = figure;
hold on
plot(V, Prec, '.', V, Psct, '.', V, Psct2D, '.',...
     ...% V, Ptrp, '.', V, Ptrp2D, '.', V, Prec + Psct + Ptrp, '.',...
     V, Prec + Psct2D + Ptrp2D, '.',...
     'MarkerSize', 20, 'LineWidth', 2.5)
xlabel('Normalized Injection Bias (eV/\Delta)', 'FontSize', 14)
ylabel('Fraction of Total Power', 'FontSize', 14)
legend('recombination', 'scattering', 'scattering above 2\Delta',...
    ...% 'trapping', 'trapping above 2\Delta', 'total',...
    'total above 2\Delta', 'Location', 'SouthEast')
title(['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0', ', ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'),...
     ...% ', ', 'v = ', num2str(vol, '%.2e'), ' \mu{m}^3'
     ])
axis tight
xlim([1, max(V)])
grid on
set(gca, 'box', 'on')
savePDF(h, 'SimPowerBudget.pdf')
savefig(h, 'SimPowerBudget.fig')
end