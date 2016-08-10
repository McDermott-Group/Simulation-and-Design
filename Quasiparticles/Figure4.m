scrsz = get(0, 'ScreenSize');
h = figure('Position', [.15, .0, .75, .9] * scrsz(4));

subplot(2, 1, 1)

load('SimNoTrap')
load('SimTrap')

nqp_no_tr(V_no_tr < 1) = [];
P_no_tr(V_no_tr < 1) = [];
nqp_tr(V_tr < 1) = [];
P_tr(V_tr < 1) = [];

hold on
plot(P_no_tr, nqp_no_tr, 'k^', P_sim_no_tr, nqp_sim_no_tr, 'k-',...
    'MarkerSize', 6, 'LineWidth', 3, 'MarkerFaceColor', 'k')
plot(P_tr, nqp_tr, 'r^', P_sim_tr, nqp_sim_tr, 'r-',...
    'MarkerSize', 6, 'LineWidth', 3, 'MarkerFaceColor', 'r')
hold off
set(gca, 'xscale', 'Log', 'yscale', 'Log')
xlabel('Power (W)', 'FontSize', 12)
ylabel('Quasiparticle Density (\mu{m}^{-3})', 'FontSize', 12)
legend({'no traps, experiment', 'no traps, simulation',...
    'with traps, experiment', 'with traps, simulation'},...
    'Location', 'SouthEast')
text(1e-11, 4.5e4, '(a)', 'FontSize', 15)
grid on
axis tight
set(gca, 'box', 'on')
ylim([5, 1.1e5])

subplot(2, 1, 2)

load('SimNIS24062016.mat')

[V, indices] = sort(V);
tau_p = tau_p(indices);

hold on
plot(E_p_n, tau_p_n, 'rh', 'MarkerSize', 6, 'LineWidth', 3,...
    'MarkerFaceColor', 'r')
plot(V, tau_p, 'r-', 'LineWidth', 3)
set(gca, 'xscale', 'Log', 'yscale', 'Log')
hold off
xlabel('Normalized Injection Bias eV/\Delta', 'FontSize', 12)
ylabel('Time Constant (\mu{s})', 'FontSize', 12)
text(2.39, 91, '(b)', 'FontSize', 15)
legend({'experiment', 'simulation'}, 'Location', 'SouthWest')
axis tight
grid on
set(gca, 'box', 'on')
ylim([25, 100])

savePDF(h, 'Figure4.pdf')