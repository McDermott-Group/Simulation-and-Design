function [t, e, n, f, n_qp, r_qp, P] = ...
    twoRegionTimeDomainModelOptimized(Tph, tspan, V, rqp, rph, c, vol, N)
%twoRegionTimeDomainModelOptimized Two-region, one that describes a normal
% metal isolator-superconductor junction (NIS) and the other one -
% a resonator, quasi-0D model for computing the quasiparticle density
% evolution.  
% 
% [t, e, n, f, n_qp, r_qp, P] = 
%   twoRegionTimeDomainModel(Tph, tspan, V, rqp, rph, c, vol, N)
%   computes the quasiparticle dynamics using a two-region quasi-0D model.
%   The quasiparticles are directly injected into the NIS junction.
%   At every point in time the phonon density is computed. This phonon
%   density is used to calculate the quasiparticle injection rate assuming
%   the pair-breaking mechanism.
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
%      per um^3, length(n_qp) == length(t),
%      r_qp is the total injection rate in n_{cp}/\tau_0, a single number,
%      P is the total injected power in W, a single number.

% Constants.
kB = 1.38064852e-23; % J / K
eV2J = 1.602176565e-19; % J / eV 

delta = 0.18e-3; % eV (aluminum superconducting gap)

% Convert all energy-related values to \Delta units.
delta = eV2J * delta; % J
Tph = kB * Tph / delta; % in units of \Delta

Tc = 1 / 1.764; % \delta/(k_B * T_c) = 1.764 at T = 0 K, BCS result

tau_0 = 438e-9; % sec
ncp = 4e6; % n_{cp} for aluminum is 4e6 um^-3
           % C. Wang et al. Nature Comm. 5, 5836 (2014)

% Maximum energy.
max_e = 2 * max(V);

% Assign the quasiparicle energies to the bins. Non-uniform energy
% spacing is implemented. To get a spacing that is close to a uniform
% set alpha to a very small positive value.
alpha = 6;
e = (1 + (max_e - 1) * sinh(alpha * (0:N) / N) / sinh(alpha))';
de = diff(e);
e = e(2:end);

% Short-hand notation.
rho_de = rho(e) .* de;
de = de';

[Gs_in, Gs_out] = Gscattering(e, de, Tph, Tc);

Gr = Grecombination(e, Tph, Tc);

Gtr = Gtrapping(e, Tph, Tc, c);

Rdirect = DirectInjection(e, rho_de, V, rqp);

[Omega1D_sct, indices_sct, N_Omega1D_no_n_de_sct, dR_no_N_Omega_de_sct] =...
    ScatteringInjection(e, Tc, Tph);

[Omega1D_rec, Gr_rec, dR_no_N_Omega_de_rec] =...
    RecombinationInjection(e, Tc, Tph);

[Omega1D_trp, indices_trp, N_Omega1D_no_n_de_trp, dR_no_N_Omega_de_trp,...
    de_trp] = TrapInjection(e, de, Tc, Tph, c);

% Initial condition for the quasiparticle distribution.
% The first half of vector n0 describes the NIS junction and the second -
% the resonator.
n0 = zeros(size(e));
n0 = [n0; n0];

% Solve the ODE. 
options = odeset('AbsTol', 1e-20);
[t, n] = ode15s(@(t, n) quasiparticleODE(t, n, Gs_in, Gs_out, Gr, Gtr,...
    Rdirect, e, de, rph, V,...
    Omega1D_sct, indices_sct, N_Omega1D_no_n_de_sct, dR_no_N_Omega_de_sct,...
    Omega1D_rec, Gr_rec, dR_no_N_Omega_de_rec,...
    Omega1D_trp, indices_trp, N_Omega1D_no_n_de_trp, dR_no_N_Omega_de_trp,...
    de_trp),...
    tspan, n0, options);

n = n(:, size(n, 2)/2+1:end);
% Occupation numbers.
f = n ./ (ones(length(t), 1) * rho_de');

% Non-equlibrium quasipartical density.
n_qp = 2 * ncp * sum(n, 2); % um^-3

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
coeff = delta * ncp * vol / tau_0;
P = coeff * P; % W
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
        rho(ej) .* Np(ei - ej, Tph) .* (ones(length(e), 1) * de) / Tc^3;
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

function Gtr = Gtrapping(e, Tph, Tc, c)
    % Simple model for the trapping matix.
    de = .5e-4; % in units of \Delta
    e_gap = de/2:de:1-de/2;
    [eg, ei] = meshgrid(e_gap, e);
    Gtr = (ei - eg).^2 .* Np(ei - eg, Tph) / Tc^3;
    Gtr = c * trapz(e_gap, Gtr, 2);
end

function R = DirectInjection(e, rho_de, V, r)
    % It is assmumed that the injection is happening at t < 0 and
    % the relaxation - at t > 0.
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

function [Omega1D, indices, N_Omega1D_no_n_de, dR_no_N_Omega_de] =...
    ScatteringInjection(e_inj, Tc, Tph)
    [ej, ei] = meshgrid(e_inj);
    N_Omega1D_no_n_de = (ei - ej).^2 .* rho(ej) .*...
        (1 - 1 ./ (ei .* ej)) .* Np(ei - ej, Tph) / Tc^3;
    Omega = ei - ej;
    Omega1D = Omega(:);
    indices = Omega1D <= 2 | N_Omega1D_no_n_de(:) <= 0;
    Omega1D(indices) = [];
    [e, Omega2D] = meshgrid(e_inj, Omega1D);
    % See Eq. (27) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
    dR_no_N_Omega_de = Omega2D.^2 .*...
                (e .* (Omega2D - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega2D - e).^2 - 1));
    dR_no_N_Omega_de(Omega2D <= e + 1) = 0;
end

function [Omega1D, Gr, dR_no_N_Omega_de] =...
    RecombinationInjection(e_inj, Tc, Tph)
    [ej, ei] = meshgrid(e_inj);
    Gr = (ei + ej).^2 .* (1 + 1 ./ (ei .* ej)) .*...
         Np(ei + ej, Tph) / Tc^3;
    Omega = ei + ej;
    Omega1D = Omega(:);

    [e, Omega2D] = meshgrid(e_inj, Omega1D);

    % See Eq. (27) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
    dR_no_N_Omega_de = Omega2D.^2 .*...
                (e .* (Omega2D - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega2D - e).^2 - 1));
    dR_no_N_Omega_de(Omega2D <= e + 1) = 0;
end

function [Omega1D, indices, N_Omega1D_no_n_de, dR_no_N_Omega_de, de] =...
    TrapInjection(e_inj, de_inj, Tc, Tph, c)
    N = length(e_inj);
    e_gap = linspace(1/(2 * N), 1 - 1 / (2 * N), N);
    [ej, ei] = meshgrid(e_gap, e_inj);

    N_Omega1D_no_n_de = c * (ei - ej).^2 .* Np(ei - ej, Tph) / Tc^3;

    Omega = ei - ej;
    Omega1D = Omega(:);
    indices = Omega1D <= 2 | N_Omega1D_no_n_de(:) <= 0;
    Omega1D(indices) = [];
    [e, Omega2D] = meshgrid(e_inj, Omega1D);
    % See Eq. (27) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
    dR_no_N_Omega_de = Omega2D.^2 .*...
                (e .* (Omega2D - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega2D - e).^2 - 1));
    dR_no_N_Omega_de(Omega2D <= e + 1) = 0;
    de = ones(size(de_inj)) / N;
end

function ndot = quasiparticleODE(t, n, Gs_in, Gs_out, Gr, Gtr, R_direct,...
        e, de, rph, V,...
    Omega1D_sct, indices_sct, N_Omega1D_no_n_de_sct, dR_no_N_Omega_de_sct,...
    Omega1D_rec, Gr_rec, dR_no_N_Omega_de_rec,...
    Omega1D_trp, indices_trp, N_Omega1D_no_n_de_trp, dR_no_N_Omega_de_trp,...
    de_trp)
    % n is in units n_{cp}, rates are in units n_{cp}/\tau_0.
    % It is assmumed that the injection is happening at t < 0 and
    % the relaxation - at t > 0.
    if t > 0
        R_direct = 0;
    end
    non_positive_n = n <= 0;
    n(non_positive_n) = 0;

    % The first half of vector n0 describes the NIS junction and
    % the second - the resonator.
    half = length(n) / 2;
    n_nis = n(1:half);
    n_res = n(half+1:end);

    % Injection due to scattering.
    if max(V) <= 3 || rph == 0
        Rph_sct = zeros(size(e));
    else
        N_Omega_sct = N_Omega1D_no_n_de_sct .* (n_nis * de);
        N_Omega1D = N_Omega_sct(:);
        N_Omega1D(indices_sct) = [];
        P2D = sum(N_Omega1D .* Omega1D_sct);
        if P2D <= 0
            Rph_sct = zeros(size(e));
        else
            dR = dR_no_N_Omega_de_sct .* (N_Omega1D * de);
            R = sum(dR)';
            Rph_sct = rph * R * P2D / sum(e .* R);
        end
    end

    % Injection due to recombination.
    if sum(n_nis) <= 0 || rph == 0
        Rph_rec = zeros(size(e));
    else
        [n_j, n_i] = meshgrid(n_nis);
        N_Omega = n_i .* Gr_rec .* n_j;

        N_Omega1D = N_Omega(:);
        P = sum(N_Omega1D .* Omega1D_rec);

        dR = dR_no_N_Omega_de_rec .* (N_Omega1D * de);
        R = sum(dR)';
        Rph_rec = rph * R * P / sum(e .* R);
    end

    % Injection due to trapping.
    if max(V) <= 2 || rph == 0
        Rph_trp = zeros(size(e));
    else
        N_Omega = N_Omega1D_no_n_de_trp .* (n_nis * de_trp);
        N_Omega1D = N_Omega(:);
        N_Omega1D(indices_trp) = [];
        P2D = sum(N_Omega1D .* Omega1D_trp);
        if P2D <= 0
            Rph_trp = zeros(size(e));
        else
            dR = dR_no_N_Omega_de_trp .* (N_Omega1D * de);
            R = sum(dR)';
            Rph_trp = rph * R * P2D / sum(e .* R);
        end
    end

    ndot_nis = Gs_in * n_nis - Gs_out .* n_nis - 2 * n_nis .* (Gr * n_nis) -...
        Gtr .* n_nis + R_direct;
    ndot_res = Gs_in * n_res - Gs_out .* n_res - 2 * n_res .* (Gr * n_res) -...
        Gtr .* n_res + Rph_rec + Rph_trp + Rph_sct;
    ndot = [ndot_nis; ndot_res];
    ndot(ndot < 0 & non_positive_n) = 0;
end