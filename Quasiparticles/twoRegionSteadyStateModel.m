function [t, e, n, f, n_qp, r_qp, P,...
    Prec, Psct, Psct2D, Ptrp, Ptrp2D,...
    Omega_rec, N_Omega_rec,...
    Omega_sct, N_Omega_sct,...
    Omega_trp, N_Omega_trp, f_nis] = ...
    twoRegionSteadyStateModel(Tph, tspan, V, rqp, rph, c,...
    vol, N, plot_flag)
% twoRegionSteadyStateModel Two-region, one corresponds to
% a normal metal-isolator-superconductor junction (NIS) and the other one -
% to a resonator, quasi-0D model for computing the steady-state
% quasiparticle densities.
%
% [t, e, n, f, n_qp, r_qp, P] = ...
%   twoRegionSteadyStateModel(Tph, tspan, V, rqp, rph, c,...
%   vol, N, plot_flag) computes the quasiparticle dynamics using
%   a two-point quasi-0D model. The quasiparticles are directly injected
%   into NIS junction. The equlibrium quasiparticle distribution and
%   phonon densities are computed. The phonon density is then used to
%   calculate the quasi-particle injection rate assuming the pair-breaking
%   mechanism.
%
%   The input parameters:
%      Tph is phonon temperature in K,
%      tspan should be in form of [ti, tf] in units of \tau_0 where ti is
%      the inital time and tf is the final time of the integration range
%     (it is assumed quasiparticles are injected at t < 0 and relaxation
%      happens at t > 0),
%      V is the applied voltage in \Delta,
%      rqp is the quasiparticle injection rate in n_{cp}/\tau_0,
%      rph is the fraction of the phonons that go to the resonator,
%      dimensionless,
%      c is the trapping (capture) rate in n_{cp}/\tau_0,
%      vol is the injection volume in um^3,
%      N is the number of the energy bins (50-1000 or so).
%
%   The output parameters:
%      t is time in \tau_0,
%      e are energies of the bins in \Delta, a vector,
%      n are the quasiparticle densities in n_{cp}/\tau_0,
%      size(n) == [length(t), length(e)],
%      f are the occupational numbers, size(f) == size(n),
%      n_qp is the non-equilibrium quasiparticle density in quasiparticles
%      per um^3,
%      length(n_qp) == length(t),
%      r_qp is the total injection rate in n_{cp}/\tau_0, a single number,
%      P is the total injected power in W, a single number.

% Constants.
kB = 1.38064852e-23; % J / K
eV2J = 1.602176565e-19; % J / eV 

delta = 0.18e-3; % eV (aluminum superconducting gap)

% Convert all energy-related values to \Delta units.
delta = eV2J * delta; % J
Tph = kB * Tph / delta; % in units of \Delta

Tc = 1 / 1.764; % \delta/(K_B * T_c) = 1.764 at T = 0 K BCS result

tau_0 = 438e-9; % sec
ncp = 4e6; % n_{cp} for aluminum is 4e6 \micro m^-3
           % C. Wang et al. Nature Comm. 5, 5836 (2014)

% Maximum energy.
max_e = 2 * max(V);

% Assign the quasiparicle energies to the bins. Non-uniform energy
% spacing is implemented. To get a spacing that is close to a uniform
% set alpha to a very small positive value.
alpha = 4;
e1 = (1 + (max_e - 1) * sinh(alpha * (0:N) / N) / sinh(alpha))';
e2 = (1 + (max_e - 1) * sinh(alpha * (1:N+1) / N) / sinh(alpha))';
de = e2 - e1;
e = (e2 + e1) / 2;

% Short-hand notation.
rho_de = rho(e) .* de;

% Bins indices the quasiparticles are injected to.
if length(V) == 2
    indices = V(1) < e & e < V(2);
else
    indices = e < V;
end

% Total injection rate.
r_qp = rqp * sum(rho_de(indices)); % n_{cp}/\tau_0

% Injection power.
P = rqp * sum(e(indices) .* rho_de(indices));
power_calib = delta * ncp * vol / tau_0;
P = power_calib * P; % W

[Gs_in, Gs_out] = Gscattering(e, de, Tph, Tc);

Gr = Grecombination(e, Tph, Tc);

Gtr = Gtrapping(e, Tph, Tc, c);

Rqp = DirectInjection(e, rho_de, V, rqp);

% Initial condition for the quasiparticle distribution.
n0 = zeros(size(e));

% Solve the ODE at the NIS junction.
options = odeset('AbsTol', 1e-10);
[t, n] = ode15s(@(t, n) quasiparticleODE(t, n,...
    Gs_in, Gs_out, Gr, Gtr, Rqp), tspan, n0, options);

n_inj = n(end, :)';

f_nis = n ./ (ones(length(t), 1) * rho_de');

% Occupation numbers.
f_inj = n_inj ./ rho_de;
f_inj(f_inj < 0) = 0;
f_inj = @(e_inj) interp1(e, f_inj, e_inj);

[Rph_rec_brute, Omega_rec, N_Omega_rec] = RecombinationInjection(e, de, f_inj,...
    V, rph, Tc, Tph);
[Rph_sct_brute, Omega_sct, N_Omega_sct] = ScatteringInjection(e, de, f_inj,...
    V, rph, Tc, Tph);
[Rph_trp_brute, Omega_trp, N_Omega_trp] = TrapInjection(e, de, f_inj,...
    V, c * rph, Tc, Tph);

Rph_rec = FastRecombinationInjection(e, de, n_inj, rph, Tc, Tph);
Rph_sct = FastScatteringInjection(e, de, n_inj, rph, Tc, Tph, V);
Rph_trp = FastTrapInjection(e, de, n_inj, rph, Tc, Tph, V, c);

n_eq  = n(end, :)';
[~, P_rec] = FastRecombinationInjection(e, de, n_eq, rph, Tc, Tph);
[~, P_sct, P_sct2D] = FastScatteringInjection(e, de, n_eq, rph, Tc, Tph, V);
[~, P_trp, P_trp2D] = FastTrapInjection(e, de, n_eq, rph, Tc, Tph, V, c);

Prec = power_calib * P_rec / P;
Psct = power_calib * P_sct / P;
Psct2D = power_calib * P_sct2D / P;
Ptrp = power_calib * P_trp / P;
Ptrp2D = power_calib * P_trp2D / P;

if plot_flag
    figure
    plot(Omega_rec, N_Omega_rec, 'LineWidth', 2)
    xlabel('Phonon Energy \Omega (\Delta)', 'FontSize', 14)
    ylabel('N_{gen}(\Omega) (n_{cp} / \Delta)', 'FontSize', 14)
    title('Phonon Occupation Numbers (Recombination)')
    set(gca, 'yscale', 'Log')
    axis tight
    grid on
    
    figure
    plot(Omega_sct, N_Omega_sct, 'LineWidth', 2)
    xlabel('Phonon Energy \Omega (\Delta)', 'FontSize', 14)
    ylabel('N_{gen}(\Omega) (n_{cp} / \Delta)', 'FontSize', 14)
    title('Phonon Occupation Numbers (Scattering)')
    set(gca, 'yscale', 'Log')
    axis tight
    grid on
    
    figure
    plot(Omega_trp, N_Omega_trp, 'LineWidth', 2)
    xlabel('Phonon Energy \Omega (\Delta)', 'FontSize', 14)
    ylabel('N_{gen}(\Omega) (n_{cp} / \Delta)', 'FontSize', 14)
    title('Phonon Occupation Numbers (Trapping)')
    set(gca, 'yscale', 'Log')
    axis tight
    grid on

    % Quaiparticle occupation numbers.
    f = n ./ (ones(length(t), 1) * rho_de');

    Pt_sct = NaN(size(n, 1), 1);
    Pt_sct2D = NaN(size(n, 1), 1);
    Pt_rec = NaN(size(n, 1), 1);
    Pt_trp = NaN(size(n, 1), 1);
    Pt_trp2D = NaN(size(n, 1), 1);
    Pt_rec_fast = NaN(size(n, 1), 1);
    Pt_sct_fast = NaN(size(n, 1), 1);
    Pt_sct2D_fast = NaN(size(n, 1), 1);
    Pt_trp_fast = NaN(size(n, 1), 1);
    Pt_trp2D_fast = NaN(size(n, 1), 1);
    power_calib = delta * ncp * vol / tau_0;
    for k = 1:size(n, 1)
        f_inj = f(k, :);
        f_inj = @(e_inj) interp1(e, f_inj, e_inj);
        [~, Omega_rec, N_Omega_rec] = RecombinationInjection(e, de,...
            f_inj, V, rph, Tc, Tph);
        [~, Omega_sct, N_Omega_sct] = ScatteringInjection(e, de,...
            f_inj, V, rph, Tc, Tph);
        [~, Omega_trp, N_Omega_trp] = TrapInjection(e, de,...
            f_inj, V, c * rph, Tc, Tph);
        Pt_sct(k) = power_calib * trapz(Omega_sct, N_Omega_sct .* Omega_sct);
        N_Omega_sct = N_Omega_sct(Omega_sct > 2);
        Omega_sct = Omega_sct(Omega_sct > 2);
        Pt_sct2D(k) = power_calib * trapz(Omega_sct, N_Omega_sct .* Omega_sct);

        Pt_rec(k) = power_calib * trapz(Omega_rec, N_Omega_rec .* Omega_rec);

        Pt_trp(k) = power_calib * trapz(Omega_trp, N_Omega_trp .* Omega_trp);
        N_Omega_trp = N_Omega_trp(Omega_trp > 2);
        Omega_trp = Omega_trp(Omega_trp > 2);
        Pt_trp2D(k) = power_calib * trapz(Omega_trp, N_Omega_trp .* Omega_trp);
    end

    figure
    Pres = rph * P;
    plot(t, Pt_sct ./ Pres, t, Pt_rec ./ Pres, t, Pt_trp ./ Pres,...
         t, (Pt_sct + Pt_rec + Pt_trp) ./ Pres,...
         t, Pt_sct2D ./ Pres, t, Pt_trp2D ./ Pres, t,...
         (Pt_sct2D + Pt_rec + Pt_trp2D) ./ Pres, 'LineWidth', 2)
    xlabel('Time (\tau_0)', 'FontSize', 14)
    ylabel('Normalized Power (P / P_{total})', 'FontSize', 14)
    legend('scattering', 'recombination', 'trapping',...
           'scattering + recombination + trapping',...
           'scattering above 2\Delta', 'trapping above 2\Delta',...
           '(scattering + recombination + trapping) above 2\Delta')
    title(['V = ',  num2str(V), ' \Delta/e (brute-force integration)']) 
    axis tight
    grid on

    for k = 1:size(n, 1)
        nt = n(k, :);
        P_sct_in = power_calib * sum(e .* (Gs_in * nt'));
        P_sct_out = power_calib * sum(e .* Gs_out .* nt');
        Pt_sct(k) = P_sct_out - P_sct_in;
        Pt_rec(k) = power_calib * sum(2 * e .* nt' .* (Gr * nt'));
        Pt_trp(k) = power_calib * sum(e .* Gtr .* nt');
    end
    
    figure
    plot(t, Pt_sct ./ P, t, Pt_rec ./ P, t, Pt_trp ./P,...
         t, (Pt_sct + Pt_rec + Pt_trp) ./ P, 'LineWidth', 2)
    xlabel('Time (\tau_0)', 'FontSize', 14)
    ylabel('Normalized Power (P / P_{total})', 'FontSize', 14)
    legend('scattering', 'recombination', 'trapping',...
           'scattering + recombination + trapping')
    title(['V = ',  num2str(V), ' \Delta/e (direct integration)']) 
    axis tight
    grid on
    
    for k = 1:size(n, 1)     
        nt = n(k, :)';
        [~, Pt_rec_fast(k)] = FastRecombinationInjection(e, de, nt, rph, Tc, Tph);
        [~, Pt_sct_fast(k), Pt_sct2D_fast(k)] = FastScatteringInjection(e, de, nt, rph, Tc, Tph, V);
        [~, Pt_trp_fast(k), Pt_trp2D_fast(k)] = FastTrapInjection(e, de, nt, rph, Tc, Tph, V, c);
    end
    
    figure
    Pt_rec_fast = power_calib * Pt_rec_fast;
    Pt_sct_fast = power_calib * Pt_sct_fast;
    Pt_trp_fast = power_calib * Pt_trp_fast;
    Pt_sct2D_fast = power_calib * Pt_sct2D_fast;
    Pt_trp2D_fast = power_calib * Pt_trp2D_fast;
    plot(t, Pt_sct_fast ./ P, t, Pt_rec_fast ./ P, t, Pt_trp_fast ./ P,...
         t, (Pt_sct_fast + Pt_rec_fast + Pt_trp_fast) ./ P,...
         t, Pt_sct2D_fast ./ P, t, Pt_trp2D_fast ./ P, t,...
         (Pt_sct2D_fast + Pt_rec_fast + Pt_trp2D_fast) ./ P, 'LineWidth', 2)
    xlabel('Time (\tau_0)', 'FontSize', 14)
    ylabel('Normalized Power (P / P_{total})', 'FontSize', 14)
    legend('scattering', 'recombination', 'trapping',...
           'scattering + recombination + trapping',...
           'scattering above 2\Delta', 'trapping above 2\Delta',...
           '(scattering + recombination + trapping) above 2\Delta')
    title(['V = ',  num2str(V), ' \Delta/e (improved calculation)']) 
    axis tight
    grid on
end

% Solve the ODE at the resonator.
options = odeset('AbsTol', 1e-10, 'RelTol', 1e-6);
[t, n] = ode15s(@(t, n) quasiparticleODE(t, n,...
    Gs_in, Gs_out, Gr, Gtr, Rph_sct + Rph_rec + Rph_trp),...
    tspan, n0, options);

if plot_flag
    figure
    plot(e, Rph_sct, e, Rph_rec, e, Rph_trp,...
        e, Rph_sct_brute, e, Rph_rec_brute, e, Rph_trp_brute, 'LineWidth', 2)
    xlabel('Energy \epsilon (\Delta)', 'FontSize', 14)
    ylabel('Injection Rate per Unit Energy (1/\tau_0\Delta)',...
        'FontSize', 14)
    legend({'R_{ph. sct.}', 'R_{ph. rec.}', 'R_{ph. trp.}',...
        'R_{ph. sct.} (brute-force)', 'R_{ph. rec.} (brute-force)',...
        'R_{ph. trp.} (brute-force)'})
    set(gca, 'yscale', 'Log')
    axis tight
    xlim([1, max(V)])
    grid on
end

% Occupational numbers.
f = n ./ (ones(length(t), 1) * rho_de');

% Non-equlibrium quasipartical density.
n_qp = 2 * ncp * sum(n, 2);
end

function normalized_density = rho(e)
    normalized_density = e ./ sqrt(e.^2 - 1);
    normalized_density(e <= 1) = 0;
end

function abs_phonon_density = Np(e, Tph)
    abs_phonon_density = 1 ./ abs(exp(-e / Tph) - 1);
end

function [Gs_in, Gs_out] = Gscattering(e, de, Tph, Tc)
% Scattering matrix defined by Eq. (6) in J. M. Martinis et al.,
% Phys. Rev. Lett. 103, 097002 (2009), as well as Eqs. (C3) and (C4) in 
% M. Lenander et al., Phys. Rev. B 84, 024501 (2011).
% These equations come from Eq. (8) in S. B. Kaplan et al.,
% Phys. Rev. B 14, 4854 (1976).
    [ej, ei] = meshgrid(e);

    Gs = (ei - ej).^2 .*...
        (1 - 1 ./ (ei .* ej)) .*...
        rho(ej) .* Np(ei - ej, Tph) .* (ones(size(e)) * de') / Tc^3;
    Gs(~isfinite(Gs)) = 0;

    Gs_in = Gs';
    Gs_out = sum(Gs, 2);
end

function Gr = Grecombination(e, Tph, Tc)
% Recombination matrix defined by Eq. (7) in J. M. Martinis et al.,
% Phys. Rev. Lett. 103, 097002 (2009), as well as Eqs. (C5) and (C6) in 
% M. Lenander et al., Phys. Rev. B 84, 024501 (2011).
% These equations come from Eq. (8) in S. B. Kaplan et al.,
% Phys. Rev. B 14, 4854 (1976).
    [ej, ei] = meshgrid(e);

    Gr = (ei + ej).^2 .*...
        (1 + 1 ./ (ei .* ej)) .*...
        Np(ei + ej, Tph) / Tc^3;
end

function Gtr = Gtrapping(e_inj, Tph, Tc, c)
    % Simple model for the trapping matix.
    N = 5 * length(e_inj);
    e_gap = linspace(0, 1, N + 1);
    [eg, ei] = meshgrid(e_gap, e_inj);
    Gtr = (ei - eg).^2 .* Np(ei - eg, Tph) / Tc^3;
    Gtr = c * trapz(e_gap, Gtr, 2);
end

function R = DirectInjection(e, rho_de, V, r)
    % It is assmumed that the injection is happening at t < 0 and
    % the relaxation at t > 0.
    % Injection voltage.
    R = zeros(size(e));
    if length(V) == 2
        % Injection into an energy band.
        R(V(1) < e & e <= V(2)) = r;
    else
        % Realistic injection.
        R(e <= V) = r;
    end
    R = R .* rho_de;
end

function [R, Omega1D, N_Omega] = ScatteringInjection(e_inj, de_inj,...
        f_inj, V, r, Tc, Tph)
    if max(V) <= 3
        R = zeros(size(e_inj));
        Omega1D = [0, 10];
        N_Omega = [0, 0];
        return
    end
    N = length(e_inj);
    Omega = linspace(0, max(V) - 1, N);
    e_initial = linspace(min(e_inj), max(e_inj), N)';
    [Omega, e] = meshgrid(Omega, e_initial);
    % See Eq. (8) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
    dN_Omega = Omega.^2 .* rho(e - Omega) .* rho(e) .* f_inj(e) .*...
        (1 - 1 ./ (e .* (e - Omega))) .* Np(Omega, Tph);
    dN_Omega(dN_Omega < 0 | ~isfinite(dN_Omega) | e - Omega < 1.01) = 0;
    N_Omega = r * trapz(e_initial, dN_Omega) / Tc^3;
    Omega1D = Omega(1, :);
    % See Eq. (27) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
    dR = Omega.^2 .* (ones(size(e_initial)) * N_Omega) .*...
                (e .* (Omega - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega - e).^2 - 1));
    dR(Omega <= e + 1) = 0;
    R = trapz(Omega(1, :), dR, 2) / Tc^3;
    R = interp1(e_initial, R, e_inj) .* de_inj;
end

function [R, Omega1D, N_Omega] = RecombinationInjection(e_inj, de_inj,...
        f_inj, V, r, Tc, Tph)
    N = 5 * length(e_inj);
    Omega = linspace(2, 2 * max(V), N);
    e_final = linspace(min(e_inj), max(e_inj), N)';
    [Omega, e] = meshgrid(Omega, e_final);
    
    % See Eq. (8) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
    dN_Omega = Omega.^2 .* rho(Omega - e) .* rho(e) .*...
        f_inj(Omega - e) .* f_inj(e) .*...
        (1 + 1 ./ (e .* (Omega - e))) .* Np(Omega, Tph);
    dN_Omega(dN_Omega < 0 | ~isfinite(dN_Omega)) = 0;
    N_Omega = r * trapz(e_final, dN_Omega) / Tc^3;

    % See Eq. (27) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
    dR = Omega.^2 .* (ones(size(e_final)) * N_Omega) .*...
                (e .* (Omega - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega - e).^2 - 1));
    dR(Omega <= e + 1) = 0;
    R = trapz(Omega(1, :), dR, 2) / Tc^3;
    R = interp1(e_final, R, e_inj) .* de_inj;
    Omega1D = Omega(1, :);
end

function [R, P, P2D] = FastScatteringInjection(e_inj, de_inj,...
        n_inj, r, Tc, Tph, V)    
    [ej, ei] = meshgrid(e_inj);
    N_Omega1D = (ei - ej).^2 .* rho(ej) .* (n_inj * de_inj') .*...
        (1 - 1 ./ (ei .* ej)) .* Np(ei - ej, Tph) / Tc^3;

    Omega = ei - ej;
    Omega1D = Omega(:);
    N_Omega1D = N_Omega1D(:);
    indices = Omega1D > 0;
    P = sum(N_Omega1D(indices) .* Omega1D(indices));
    indices = Omega1D <= 2 | N_Omega1D <= 0;
    N_Omega1D(indices) = [];
    Omega1D(indices) = [];
    P2D = sum(N_Omega1D .* Omega1D);
    
    if max(V) <= 3 || r == 0 || isempty(Omega1D)
        R = zeros(size(e_inj));
        return
    end

    [e, Omega2D] = meshgrid(e_inj, Omega1D);

    % See Eq. (27) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
    dR = Omega2D.^2 .* (N_Omega1D * de_inj') .*...
                (e .* (Omega2D - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega2D - e).^2 - 1));
    dR(Omega2D <= e + 1) = 0;
    R = sum(dR)';
    R = r * R * P2D / sum(e_inj .* R);
end

function [R, P] = FastRecombinationInjection(e_inj, de_inj,...
        n_inj, r, Tc, Tph)
    [ej, ei] = meshgrid(e_inj);
    
    Gr = (ei + ej).^2 .* (1 + 1 ./ (ei .* ej)) .*...
        Np(ei + ej, Tph) / Tc^3;
    n_inj(n_inj < 0) = 0;
    [n_j, n_i] = meshgrid(n_inj);
    N_Omega = n_i .* Gr .* n_j;
    
    Omega = ei + ej;
    Omega1D = Omega(:);
    N_Omega1D = N_Omega(:);
    P = sum(N_Omega(:) .* Omega1D);

    [e, Omega2D] = meshgrid(e_inj, Omega1D);

    % See Eq. (27) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
    dR = Omega2D.^2 .* (N_Omega1D * de_inj') .*...
                (e .* (Omega2D - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega2D - e).^2 - 1));
    dR(Omega2D <= e + 1) = 0;
    R = sum(dR)';
    R = r * R * P / sum(e_inj .* R);
end

function [R, P, P2D] = FastTrapInjection(e_inj, de_inj,...
        n_inj, r, Tc, Tph, V, c)
    N = length(e_inj);
    e_gap = linspace(1/(2 * N), 1 - 1 / (2 * N), N);
    
    [ej, ei] = meshgrid(e_gap, e_inj);
    N_Omega = c * (ei - ej).^2 .* (n_inj * ones(size(e_gap)) / N) .*...
                Np(ei - ej, Tph) / Tc^3;
            
    Omega = ei - ej;
 
    Omega1D = Omega(:);
    N_Omega1D = N_Omega(:);
    P = sum(N_Omega1D .* ei(:));
    indices = Omega1D <= 2 | N_Omega1D <= 0;
    N_Omega1D(indices) = [];
    Omega1D(indices) = [];
    P2D = sum(N_Omega1D .* Omega1D);
    
    if max(V) <= 2 || r == 0 || c == 0 || isempty(Omega1D)
        R = zeros(size(e_inj));
        return
    end

    [e, Omega2D] = meshgrid(e_inj, Omega1D);

    % See Eq. (27) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
    dR = Omega2D.^2 .* (N_Omega1D * de_inj') .*...
                (e .* (Omega2D - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega2D - e).^2 - 1));
    dR(Omega2D <= e + 1) = 0;
    R = sum(dR)';
    R = r * R * P2D / sum(e_inj .* R);
end

function [R, Omega1D, N_Omega] = TrapInjection(e_inj, de_inj,...
        f_inj, V, r, Tc, Tph)
    if max(V) <= 2
        R = zeros(size(e_inj));
        Omega1D = [0, 10];
        N_Omega = [0, 0];
        return
    end
    N = 5 * length(e_inj);
    Omega = linspace(0, 2 * max(V), N);
    e_inital = linspace(min(e_inj), max(e_inj), N)';
    [Omega, e] = meshgrid(Omega, e_inital);
    
    % See Eq. (8) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976),
    % assuming \Delta=0.
    dN_Omega = Omega.^2 .* rho(e) .* f_inj(e) .* Np(Omega, Tph);
    dN_Omega(dN_Omega < 0 | ~isfinite(dN_Omega) |...
        e - Omega >= 1 | e - Omega < 0) = 0;
    N_Omega = r * trapz(e_inital, dN_Omega) / Tc^3;
    
    % See Eq. (27) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
    dR = Omega.^2 .* (ones(size(e_inital)) * N_Omega) .*...
                (e .* (Omega - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega - e).^2 - 1));
    dR(Omega <= e + 1) = 0;
    R = trapz(Omega(1, :), dR, 2) / Tc^3;
    R = interp1(e_inital, R, e_inj) .* de_inj;
    Omega1D = Omega(1, :);
end

function ndot = quasiparticleODE(t, n, Gs_in, Gs_out, Gr, Gtr, R)
    % n is in units n_{cp}, rates are in units n_{cp}/\tau_0.
    % It is assmumed that the injection is happening at t < 0 and
    % the relaxation - at t > 0.
    if t > 0
        R = 0;
    end
    non_positive_n = n <= 0;
    ndot = Gs_in * n - Gs_out .* n - 2 * n .* (Gr * n) - Gtr .* n + R;
    ndot(ndot < 0 & non_positive_n) = 0;
end