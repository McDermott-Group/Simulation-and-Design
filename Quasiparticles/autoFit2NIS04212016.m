function autoFit2NIS04212016
%autoFit2NIS04212016 Fitting to the `NIS04212016` dataset.

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 2.232e-05; r_phonon = 5.003e-01; c = 5.405e-02; vol = 2.600e+04;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 250;

% Time domain data fit.
data = load('NIS04212016.mat');
E_p_n = data.NearTrapPoisoning(:, 2);
nqp_p_n = data.NearTrapPoisoning(:, 4);

E_r_n = data.NearTrapRecovery(:, 2);
nqp_r_n = data.NearTrapRecovery(:, 4);

V = [E_p_n; E_r_n];
nqp = [nqp_p_n; nqp_r_n];

options = optimset('Display', 'iter', 'MaxIter', 100,...
    'TolFun', 1e-2);
x = fminsearch(@(x) simulations(x, Tph, tspan, V, nqp, vol, N),...
    [r_direct, r_phonon, c], options);

disp(['r_direct = ', num2str(x(1), '%.3e'), '; ',...
      'r_phonon = ', num2str(x(2), '%.3e'), '; ',...
      'c = ', num2str(x(3), '%.3e'), '; ',...
      'vol = ', num2str(vol, '%.3e'), ';'])
end

function error = simulations(x, Tph, tspan, V, nqp, vol, N)
    indices = (V > 1) & (nqp > 0);
    nqp = nqp(indices);
    V = V(indices);
    nqp_sim = NaN(size(V));
    r_qp = x(1);
    r_ph = x(2);
    c = x(3);
    parfor k = 1:length(V)
        [~, ~, ~, ~, n_qp] = twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V(k), r_qp, r_ph, c, vol, N);
        nqp_sim(k) = n_qp(end);
    end
    error = sum((log(nqp_sim ./ nqp)).^2);
end