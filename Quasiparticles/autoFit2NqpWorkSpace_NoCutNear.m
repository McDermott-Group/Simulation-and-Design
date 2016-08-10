function autoFit2NqpWorkSpace_NoCutNear
%autoFit2NqpWorkSpace_NoCutNear Fitting to the `NqpWorkSpace` dataset.

% With traps.
% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 5.137e-05; r_phonon = 1.462e+00; c = 1.738e-02; vol = 1.669e+04;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 200;

delta = 0.18e-3; % eV (aluminum superconducting gap)
data = load('NqpWorkSpace.mat');

% With traps.
V = data.NoCutNear(:, 1) / delta;
P = data.NoCutNear(:, 2);
nqp = data.NoCutNear(:, 4);

options = optimset('Display', 'iter', 'MaxIter', 30, 'TolFun', 1e-2);
x = fminsearch(@(x) simulations(x, Tph, tspan, V, P, nqp, N),...
    [r_direct, r_phonon, c, vol], options);

disp(['r_direct = ', num2str(x(1), '%.3e'), '; ',...
      'r_phonon = ', num2str(x(2), '%.3e'), '; ',...
      'c = ', num2str(x(3), '%.3e'), '; ',...
      'vol = ', num2str(x(4), '%.3e'), ';'])
end

function error = simulations(x, Tph, tspan, V, P, nqp, N)
    indices = (V > 1.1) & (nqp > 0) & (V < 4);
    P = P(indices);
    nqp = nqp(indices);
    V = V(indices);
    nqp_sim = NaN(size(V));
    P_sim = NaN(size(V));
    r_qp = x(1);
    r_ph = x(2);
    c = x(3);
    vol = x(4);
    parfor k = 1:length(V)
        [~, ~, ~, ~, n_qp, ~, P_sim(k)] = ...
            twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V(k), r_qp, r_ph, c, c, vol, N);
        nqp_sim(k) = n_qp(end);
    end
    error = sum((log(nqp_sim ./ nqp)).^2 +...
                (log(P_sim ./ P)).^2);
end