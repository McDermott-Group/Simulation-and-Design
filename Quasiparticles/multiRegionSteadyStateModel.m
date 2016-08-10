function [t, e, n, f, n_qp, r_qp, P] = ...
    multiRegionSteadyStateModel(Tph, tspan, V, rqp, rph,...
    cqp, cph, vol, N)
% multiRegionSteadyStateModel Multi-region, quasi-0D model for computing
% the steady-state quasiparticle densities.
%
% t, e, n, f, n_qp, r_qp, P,...
%     Gamma_tr_inj, Gamma_tr_ph, Gamma_r_inj, Gamma_r_ph,...
%     n_nis, f_nis, n_qp_nis] = ...
%     twoRegionSteadyStateModelOptimized(Tph, tspan, V, rqp, rph,...
%     cqp, cph, vol, N) computes the quasiparticle dynamics using
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
%      cqp is the trapping efficiency at the NIS junction (dimensionless),
%      cph is the trapping efficiency at the resonator (dimensionless),
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

Tc = 1 / 1.764; % \delta/(K_B * T_c) = 1.764 at T = 0 K BCS result

tau_0 = 438e-9; % sec
ncp = 4e6; % n_{cp} for aluminum is 4e6 um^-3
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

Gtr = Gtrapping(e, Tph, Tc, cqp);

Rqp = DirectInjection(e, rho_de, V, rqp);

% Initial condition for the quasiparticle distribution.
n0 = zeros(size(e));

for iter = 1:4
    % Solve the ODE at the NIS junction.
    options = odeset('AbsTol', 1e-10);
    [t, n] = ode15s(@(t, n) quasiparticleODE(t, n,...
        Gs_in, Gs_out, Gr, Gtr, Rqp), tspan, n0, options);

    % Equilibrium distribution.
    n_inj = n(end, :)';

    Rph_rec = RecombinationInjection(e, de, n_inj, rph, Tc, Tph);
    Rph_sct = ScatteringInjection(e, de, n_inj, rph, Tc, Tph, V);
    Rph_trp = TrapInjection(e, de, n_inj, rph, Tc, Tph, V, cph);

    Rqp = Rph_sct + Rph_rec + Rph_trp;
end

% Occupation numbers.
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

function R = ScatteringInjection(e_inj, de_inj, n_inj, r, Tc, Tph, V)
    if max(V) <= 3 || r == 0
        R = zeros(size(e_inj));
        return
    end

    [ej, ei] = meshgrid(e_inj);
    N_Omega1D = (ei - ej).^2 .* rho(ej) .* (n_inj * de_inj') .*...
        (1 - 1 ./ (ei .* ej)) .* Np(ei - ej, Tph) / Tc^3;

    Omega = ei - ej;
    Omega1D = Omega(:);
    N_Omega1D = N_Omega1D(:);
    indices = Omega1D <= 2 | N_Omega1D <= 0;
    N_Omega1D(indices) = [];
    Omega1D(indices) = [];
    P2D = sum(N_Omega1D .* Omega1D);
    
    if isempty(Omega1D)
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

function R = RecombinationInjection(e_inj, de_inj, n_inj, r, Tc, Tph)
    [ej, ei] = meshgrid(e_inj);
    
    Gr = (ei + ej).^2 .* (1 + 1 ./ (ei .* ej)) .* Np(ei + ej, Tph) / Tc^3;
    n_inj(n_inj < 0) = 0;
    [n_j, n_i] = meshgrid(n_inj);
    N_Omega = n_i .* Gr .* n_j;
    
    Omega = ei + ej;
    Omega1D = Omega(:);
    N_Omega1D = N_Omega(:);
    P = sum(N_Omega1D .* Omega1D);

    [e, Omega2D] = meshgrid(e_inj, Omega1D);

    % See Eq. (27) in S. B. Kaplan et al., Phys. Rev. B 14, 4854 (1976).
    dR = Omega2D.^2 .* (N_Omega1D * de_inj') .*...
                (e .* (Omega2D - e) + 1) ./...
                (sqrt(e.^2 - 1) .* sqrt((Omega2D - e).^2 - 1));
    dR(Omega2D <= e + 1) = 0;
    R = sum(dR)';
    R = r * R * P / sum(e_inj .* R);
end

function R = TrapInjection(e_inj, de_inj, n_inj, r, Tc, Tph, V, c)
    if max(V) <= 2 || r == 0 || c == 0
        R = zeros(size(e_inj));
        return
    end

    N = length(e_inj);
    e_gap = linspace(1/(2 * N), 1 - 1 / (2 * N), N);
    
    [ej, ei] = meshgrid(e_gap, e_inj);
    N_Omega = c * (ei - ej).^2 .* (n_inj * ones(size(e_gap)) / N) .*...
                Np(ei - ej, Tph) / Tc^3;
            
    Omega = ei - ej;
 
    Omega1D = Omega(:);
    N_Omega1D = N_Omega(:);
    indices = Omega1D <= 2 | N_Omega1D <= 0;
    N_Omega1D(indices) = [];
    Omega1D(indices) = [];
    P2D = sum(N_Omega1D .* Omega1D);
    
    if isempty(Omega1D)
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