function steadyStateModelPhononTimes
%steadyStateModelPhononTimes Explore the effective phonon lifetimes
%at different injection voltages.

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon_n and r_phonon_f dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 1.127e-04; r_phonon = 3.069e-01; c = 9.030e-02; vol = 2.600e+04;

V = linspace(1.01, 7, 25); % in units of \delta

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

Prec = nan(size(V));
Psct = nan(size(V));
Ptrp = nan(size(V));
tau_B_rec = nan(size(V));
tau_B_sct = nan(size(V));
tau_B_trp = nan(size(V));
parfor kV = 1:length(V)
    Vsim = V(kV);
    [~, e, ~, ~, ~, ~, ~,...
        Prec(kV), Psct(kV), ~, Ptrp(kV), ~,...
        Omega_rec, N_Omega_rec,...
        Omega_sct, N_Omega_sct,...
        Omega_trp, N_Omega_trp, f_nis] =...
            twoRegionSteadyStateModel(Tph, tspan,...
            Vsim, r_direct, r_phonon, c, vol, N, false);
    f_nis = f_nis(end, :);
    f_nis(f_nis < 0 | ~isfinite(f_nis)) = 0;
    tau_B_rec(kV) = tau_ph(e, f_nis, Omega_rec, N_Omega_rec);
    tau_B_sct(kV) = tau_ph(e, f_nis, Omega_sct, N_Omega_sct);
    tau_B_trp(kV) = tau_ph(e, f_nis, Omega_trp, N_Omega_trp);

    fprintf('*')
end

tau_B = (tau_B_rec .* Prec +  tau_B_sct .* Psct + tau_B_trp .* Ptrp) ./...
        (Prec + Psct + Ptrp);

h = figure;
hold on
plot(V, tau_B_rec, V, tau_B_sct, V, tau_B_trp, V, tau_B, 'LineWidth', 3)
xlabel('Normalized Injection Bias (eV/\Delta)', 'FontSize', 14)
ylabel('Effective Phonon Lifetime \tau^{ph} (ns)', 'FontSize', 14)
legend('recombination', 'scattering', 'trapping', 'weighted average',...
    'Location', 'SouthEast')
title(['r_{qp} = ', num2str(r_direct, '%.2e'), '/\tau_0', ', ',...
     'r_{ph} = ', num2str(r_phonon, '%.3f'), ', ',...
     'c = ', num2str(c, '%.3f'), ', ',...
     'v = ', num2str(vol, '%.2e'), ' \mu{m}^3'])
axis tight
grid on
set(gca, 'yscale', 'Log')
set(gca, 'box', 'on')
savePDF(h, 'SimAveragePhononTimes.pdf')
end

function tau = tau_ph(e_comp, f, Omega, N_Omega)
    % See Table II in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
    tau0_ph = .242; % ns

    min_e = min(e_comp);
    N = 10000;
    tau_B = NaN(size(Omega));
    for k = 1:length(Omega)
        if Omega(k) - 1 > min_e
            e = linspace(min_e, Omega(k) - min_e, N);
            f1 = interp1(e_comp, f, e);
            f1(~isfinite(f1)) = 0;
            f2 = interp1(e_comp, f, Omega(k) - e);
            f2(~isfinite(f2)) = 0;
            max(f1)
            max(f2)
            tau_B(k) = tau0_ph * pi / trapz(e,...
                (e .* (Omega(k) - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega(k) - e).^2 - 1)) .*...
                (1 - f1 - f2));
        else
            tau_B(k) = 0;
        end
    end
    tau = sum(tau_B .* N_Omega) / sum(N_Omega);
end