%simple simulation of a cavity trajectory
clear all;
f_c = 4.9135e9;
f_10 = 4.55e9;
chi = 2.14e6;
kappa = 6.78e6;
alpha = -300e6;

n_d = 100;
drive = sqrt(kappa^2*n_d/4);
time = 0:0.5:2000;
Nmax = 5e5;


%[chi0, chi1] = NumericalCavityTLS(f_c, f_10, chi, [], 0, 0, 5e5, 1)
[chi0, chi1, ~] = NumericalCavityMLS(f_c, f_10, alpha, chi, [], 0, 0, 5e5, 1);
f_dr = f_c;

n = 1:Nmax-1;

cavitySimulationChiShift(time, drive, f_dr, n(2:end), chi0(2:end), chi1(2:end), f_c, kappa, 1);