function autoFit2TrapNoTrapFixedVolume_NoTrap
%autoFit2TrapNoTrapFixedVolume_NoTrap Fitting to the `TrapNoTrap` dataset.

% No traps.
% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 1.501e-05; r_phonon = 1.076e+00; c = 1.159e-02; vol = 2.600e+04;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 250;

delta = 0.18e-3; % eV (aluminum superconducting gap)
data = load('TrapNoTrap.mat');

% No traps.
V = data.NoTrap(:, 5) / delta;
P = data.NoTrap(:, 6);
nqp = data.NoTrap(:, 8) - min(data.NoTrap(:, 8));

options = optimset('Display', 'iter', 'MaxIter', 100,...
    'TolFun', 1e-2);
x = fminsearch(@(x) simulations(x, Tph, tspan, V, P, nqp, vol, N),...
    [r_direct, r_phonon, c], options);

disp(['r_direct = ', num2str(x(1), '%.3e'), '; ',...
      'r_phonon = ', num2str(x(2), '%.3e'), '; ',...
      'c = ', num2str(x(3), '%.3e'), '; ',...
      'vol = ', num2str(vol, '%.3e'), '; '])
end

function error = simulations(x, Tph, tspan, V, P, nqp, vol, N)
    indices = (V > 1) & (nqp > 0) & (V < 4);
    P = P(indices);
    nqp = nqp(indices);
    V = V(indices);
    nqp_sim = NaN(size(V));
    P_sim = NaN(size(V));
    r_qp = x(1);
    r_ph = x(2);
    c = x(3);
    parfor k = 1:length(V)
        [~, ~, ~, ~, n_qp, ~, P_sim(k)] = ...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V(k), r_qp, r_ph, c, vol, N);
        nqp_sim(k) = n_qp(end);
    end
    maxqp = max(nqp);
    maxP = max(P);
    error = sum(log(((nqp_sim - nqp) / maxqp).^2) +... 
                log(((P_sim - P) / maxP).^2));
end