function autoFit2NqpWorkSpaceVariableBias_CutNear
%autoFit2NqpWorkSpaceVariableBias_CutNear Fitting to the `NqpWorkSpace` dataset.

% With traps.
% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 2.535e-05; r_phonon = 1.697e+00; c = 1.444e-02; vol = 1.968e+04;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 100;

delta = 0.184e-3; % eV (aluminum superconducting gap)
data = load('NqpWorkSpace.mat');

% With traps.
V = data.CutNear(:, 1);
P = data.CutNear(:, 2);
nqp = data.CutNear(:, 4);

options = optimset('Display', 'iter', 'MaxIter', 30, 'TolFun', 1e-2);
x = fminsearch(@(x) simulations(x, Tph, tspan, V, P, nqp, N),...
    [r_direct, r_phonon, c, vol, 0], options);

disp(['r_direct = ', num2str(x(1), '%.3e'), '; ',...
      'r_phonon = ', num2str(x(2), '%.3e'), '; ',...
      'c = ', num2str(x(3), '%.3e'), '; ',...
      'vol = ', num2str(x(4), '%.3e'), '; ',...
      'delta = ', num2str(delta + x(5)^2, '%.4e'), ';'])
end

function error = simulations(x, Tph, tspan, V, P, nqp, N)
    delta = 0.18e-3; % eV (aluminum superconducting gap)
    V = V / (delta + (x(5))^2);
    indices = (V > 1) & (nqp > 0);
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