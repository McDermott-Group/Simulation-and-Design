function steadyStateModelPowerLawExponent
%steadyStateModelPowerLawExponent Explore the quasiparticle steady-state
% density at different injection rates (resistances).

r_direct = logspace(-1, 1, 5) * 1e-05; % in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
r_phonon = .5; % dimensionless
c = [0, 0.005, .01, .02, .05, .1, .15, .2]; % dimensionless
vol = 2.6e+04; % um^3
V = linspace(1.01, 1.25, 20); % in units of \delta

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

nnis = NaN(length(V), length(r_direct));
nres = NaN(length(V), length(r_direct));
P = NaN(size(nres));
power_nis = NaN(size(r_direct));
power_res = NaN(size(r_direct));
power_nis_mean = NaN(size(c));
power_nis_std = NaN(size(c));
power_res_mean = NaN(size(c));
power_res_std = NaN(size(c));
for kc = 1:length(c)
    for krqp = 1:length(r_direct)
        rqpsim = r_direct(krqp);
        csim = c(kc);
        parfor kV = 1:length(V)
            Vsim = V(kV);
            [~, ~, ~, ~, n_res, ~, P(kV, krqp),...
             ~, ~, ~, ~, ~, ~, n_nis] = ...
                twoRegionSteadyStateModelOptimized(Tph, tspan,...
                Vsim, rqpsim, r_phonon, csim, csim, vol, N);
            nnis(kV, krqp) = max(n_nis);
            nres(kV, krqp) = max(n_res);
            fprintf('*')
        end
        pnis = polyfit(log(P(:, krqp)), log(nnis(:, krqp)), 1);
        pres = polyfit(log(P(:, krqp)), log(nres(:, krqp)), 1);
        power_nis(krqp) = pnis(1);
        power_res(krqp) = pres(1);
    end
    power_nis_mean(kc) = mean(power_nis);
    power_nis_std(kc) = std(power_nis);
    power_res_mean(kc) = mean(power_res);
    power_res_std(kc) = std(power_res);
end

h = figure;
hold on
errorbar(c, power_nis_mean, power_nis_std, '.',...
         'MarkerSize', 10, 'LineWidth', 2)
errorbar(c, power_res_mean, power_res_std, '.',...
         'MarkerSize', 10, 'LineWidth', 2)
hold off
xlabel('Trapping Strength c', 'FontSize', 14)
ylabel('Power-Law Exponent \alpha', 'FontSize', 14)
title('n_{qp} \propto P^{\alpha}')
legend('NIS junction', 'resonator', 'Location', 'SouthEast')
axis tight
ylim([0 2])
grid on
set(gca, 'box', 'on')
savePDF(h, 'SimPowerLawExponent.pdf')
end