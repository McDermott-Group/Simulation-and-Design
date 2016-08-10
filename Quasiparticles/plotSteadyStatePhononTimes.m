function plotSteadyStatePhononTimes
%plotSteadyStatePhononTimes Phonon lifetimes.

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon_n and r_phonon_f dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 1.127e-04; r_phonon = 3.069e-01; c = 9.030e-02; vol = 2.600e+04;

N = 250;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% See Table II in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
tau0_ph = .242; % ns

V = 1.5;

[~, e_comp, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, f_nis] =...
    twoRegionSteadyStateModelOptimized(Tph, tspan, V,...
    r_direct, r_phonon, c, c, vol, N);

Omega = 0:.05:(2*V);
min_e = min(e_comp);
N = 5000;
tau_B = NaN(size(Omega));
tau_s = NaN(size(Omega));
f_nis = f_nis(end, :);
f_nis(f_nis < 0 | ~isfinite(f_nis)) = 0;
for k = 1:length(Omega)
    if Omega(k) - 1 > min_e
        e = linspace(min_e, Omega(k) - min_e, N);
        tau_B(k) = tau0_ph * pi / trapz(e,...
            (e .* (Omega(k) - e) + 1) ./...
            (sqrt(e.^2 - 1) .* sqrt((Omega(k) - e).^2 - 1)) .*...
            (1 - interp1(e_comp, f_nis, e) -...
                 interp1(e_comp, f_nis, Omega(k) - e)));
    else
        tau_B(k) = 0;
    end
    e = linspace(min_e, 2*V, N);
    f1 = interp1(e_comp, f_nis, e);
    f1(~isfinite(f1)) = 0;
    f2 = interp1(e_comp, f_nis, Omega(k) + e);
    f2(~isfinite(f2)) = 0;
    tau_s(k) = tau0_ph * (pi / 2) / trapz(e,...
        (e .* (Omega(k) + e) - 1) ./...
        (sqrt(e.^2 - 1) .* sqrt((Omega(k) + e).^2 - 1)) .*...
        (f1 - f2));
end
h = figure;
semilogy(Omega, tau_B, Omega, tau_s, 'LineWidth', 3)
xlabel('Phonon Energy \Omega/\Delta', 'FontSize', 14)
ylabel('Phonon Lifetime \tau^{ph} (ns)', 'FontSize', 14)
legend({'pair-breaking \tau^{ph}_B', 'scattering \tau^{ph}_s'},...
    'Location', 'SouthWest')
axis tight
grid on

savePDF(h, 'SimPhononTimes.pdf')

h = figure;
semilogy(Omega, tau_B, 'LineWidth', 3)
xlabel('Phonon Energy \Omega/\Delta', 'FontSize', 14)
ylabel('Phonon Lifetime \tau^{ph}_B (ns)', 'FontSize', 14)
title('Pair-Breaking \tau^{ph}_B')
axis tight
grid on

savePDF(h, 'SimPhononTimesPairBreaking.pdf')

end