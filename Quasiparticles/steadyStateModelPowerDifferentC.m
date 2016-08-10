function steadyStateModelPowerDifferentC
%steadyStateModelPowerDifferentC Explore the quasiparticle power budget
%at different trapping rates.

r_direct = 1e-05; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon = 1; % dimensionless
c = linspace(0, .2, 21); % dimensionless
vol = 2.6e+04; % um^3
V = 1.25; % in units of \delta

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

Prec = nan(size(c));
Psct = nan(size(c));
Psct2D = nan(size(c));
Ptrp = nan(size(c));
Ptrp2D = nan(size(c));
parfor kc = 1:length(c)
    csim = c(kc);
    [~, ~, ~, ~, ~, ~, ~, ...
        Prec(kc), Psct(kc), Psct2D(kc), Ptrp(kc), Ptrp2D(kc)] =...
        twoRegionSteadyStateModel(Tph, tspan,...
        V, r_direct, r_phonon, csim, vol, N, false);
    fprintf('*')
end

h = figure;
hold on
plot(c, Prec, c, Psct, c, Psct2D, c, Ptrp, c, Ptrp2D,...
    c, Prec + Psct + Ptrp, c, Prec + Psct2D + Ptrp2D,...
    'MarkerSize', 10, 'LineWidth', 2)
xlabel('Trapping Strength c', 'FontSize', 14)
ylabel('Fraction of Total Power', 'FontSize', 14)
legend('recombination', 'scattering', 'scattering above 2\Delta',...
    'trapping', 'trapping above 2\Delta', 'total',...
    'total above 2\Delta', 'Location', 'SouthEast')
title(['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0', ', ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'v = ', num2str(vol, '%.2e'), ' \mu{m}^3, ',...
     'eV/\Delta = ', num2str(V, '%.3f'), 'V'])
axis tight
grid on
set(gca, 'box', 'on')
savePDF(h, 'SimPowerBudget.pdf')
end