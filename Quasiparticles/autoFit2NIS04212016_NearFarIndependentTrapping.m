function autoFit2NIS04212016_NearFarIndependentTrapping
%autoFit2NIS04212016_NearFarIndependentTrapping Fitting to
% the `NIS04212016` dataset.

% r_direct in units of 1/\tau_0, assuming n_{qp} in units of n_{cp}
% r_phonon_n and r_phonon_f dimensionless
% c dimensionless
% vol in units of um^3
r_direct = 1.679e-04; r_phonon_n = 3.044e-01; r_phonon_f = 7.463e-02; c_qp = 7.463e-02; c_n = 1.705e-01; c_f = 3.499e-02; vol = 2.600e+04;

Tph = 0.051; % K
tspan = [-510, -10]; % in units of \tau_0

% Number of the energy bins.
N = 250;

% Data.
data = load('NIS04212016.mat');
E_p_n = data.NearTrapPoisoning(:, 2);
nqp_p_n = data.NearTrapPoisoning(:, 4);

E_r_n = data.NearTrapRecovery(:, 2);
nqp_r_n = data.NearTrapRecovery(:, 4);

V_n = [E_p_n; E_r_n];
nqp_n = [nqp_p_n; nqp_r_n];

E_p_f = data.FarTrapPoisoning(:, 2);
nqp_p_f = data.FarTrapPoisoning(:, 4);

E_r_f = data.FarTrapRecovery(:, 2);
nqp_r_f = data.FarTrapRecovery(:, 4);

V_f = [E_p_f; E_r_f];
nqp_f = [nqp_p_f; nqp_r_f];

options = optimset('Display', 'iter', 'MaxIter', 125,...
    'TolFun', 1e-2);
x = fminsearch(@(x) simulations(x, Tph, tspan,...
    V_n, nqp_n, V_f, nqp_f, vol, N),...
    [r_direct, r_phonon_n, r_phonon_f, c_qp, c_n, c_f], options);

disp(['r_direct = ', num2str(x(1), '%.3e'), '; ',...
      'r_phonon_n = ', num2str(x(2), '%.3e'), '; ',...
      'r_phonon_f = ', num2str(x(3), '%.3e'), '; ',...
      'c_qp = ', num2str(x(3), '%.3e'), '; ',...
      'c_n = ', num2str(x(4), '%.3e'), '; ',...
      'c_f = ', num2str(x(5), '%.3e'), '; ',...
      'vol = ', num2str(vol, '%.3e'), ';'])
end

function error = simulations(x, Tph, tspan, V_n, nqp_n, V_f, nqp_f, vol, N)
    r_qp = x(1);
    r_ph_n = x(2);
    r_ph_f = x(3);
    c_qp = x(4);
    c_n = x(5);
    c_f = x(6);

    indices = (V_n > 1) & (nqp_n > 0);
    nqp_n = nqp_n(indices);
    V_n = V_n(indices);
    nqp_sim_n = NaN(size(V_n));

    parfor k = 1:length(V_n)
        [~, ~, ~, ~, n_qp] = twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V_n(k), r_qp, r_ph_n, c_qp, c_n, vol, N);
        nqp_sim_n(k) = n_qp(end);
    end
    
    indices = (V_f > 1) & (nqp_f > 0);
    nqp_f = nqp_f(indices);
    V_f = V_f(indices);
    nqp_sim_f = NaN(size(V_f));
    parfor k = 1:length(V_f)
        [~, ~, ~, ~, n_qp] = twoRegionSteadyStateModelOptimized(Tph, tspan,...
            V_f(k), r_qp, r_ph_f, c_qp, c_f, vol, N);
        nqp_sim_f(k) = n_qp(end);
    end
    error = sum((log(nqp_sim_f ./ nqp_f)).^2 +...
                (log(nqp_sim_n ./ nqp_n)).^2);
end