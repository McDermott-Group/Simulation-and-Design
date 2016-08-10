function plotSingleSolution
%plotSingleSolution Quasiparticle dynamics plots.

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 2e-05; r_phonon = 5e-01; c = 5e-02; vol = 5.000e+03; 

N = 50;

Tph = 0.051; % K
tspan = [-500, 2000]; % in units of \tau_0

V = 2;

% [t, e, n, f, n_qp] = ...
%     twoRegionSteadyStateModel(Tph, tspan, V,...
%     r_direct, r_phonon, c, vol, N, true);
% [t, e, n, f, n_qp] = ...
%     twoRegionSteadyStateModelOptimized(Tph, tspan, V,...
%     r_direct, r_phonon, c, vol, N);
% clear twoRegionTimeDomainModel
% [t, e, n, f, n_qp] = ...
%     twoRegionTimeDomainModel(Tph, tspan, V,...
%     r_direct, r_phonon, c, vol, N);
[t, e, n, f, n_qp] = ...
    twoRegionTimeDomainModelOptimized(Tph, tspan, V,...
    r_direct, r_phonon, c, vol, N);

figure
plot(t, n_qp, 'LineWidth', 3)
hold on
xlabel('Time (\tau_0)', 'FontSize', 14)
ylabel('n_{\rm qp} (\mu m^{-3})', 'FontSize', 14)
title({'Quasipaticle Dynamics',...
       '(injection at t < 0, recovery at t > 0)'})
grid on
grid minor
axis tight

figure
[Ind1, Ind2] = ndgrid(t, e);
hndl = surf(Ind1, Ind2, f);
set(gca, 'View', [0 90])
set(hndl, 'LineStyle', 'none', 'FaceColor', 'interp', 'FaceLighting', 'phong');
axis tight
colormap(jet)
colorbar
xlabel('Time (\tau_0)', 'FontSize', 14)
ylabel('Energy (\epsilon/\Delta)', 'FontSize', 14)
title({'Occupational Number f(\epsilon) Time Evolution',...
       '(injection at t < 0, recovery at t > 0)'})

figure
n(n < 0) = NaN;
f(f < 0) = NaN;
semilogy(e, n(end, :), e, f(end, :), 'LineWidth', 3)
xlabel('Energy (\Delta)', 'FontSize', 14)
ylabel('n(\epsilon), f(\epsilon)', 'FontSize', 14)
legend('n(\epsilon)', 'f(\epsilon)')
axis tight
xlim([1, max(V)])
grid on

[tau_p, err_p, tau_r, err_r] = estimateTimeConstants(t, n_qp, true)

end