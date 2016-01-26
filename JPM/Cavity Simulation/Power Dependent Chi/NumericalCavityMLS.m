function [chi0, chi1, chi2] = NumericalCavityMLS(f_c, f_10, alpha, ...
                                chi_lp, epsilon, f_d, kappa, Nmax, plots)
%%NUMERICALCAVITYTLS -- Calculate chi shifts and cavity state as a function
%%of drive for a three level qubit. The cavity is treated semiclasically, and
%%the qubit is treated using an exact diagonalization of the
%%Jaynes-Cummings hamiltonian in the total excitation number basis. See
%%also NUMERICALCAVITYTLS and these additional references:
%%Larson, Phys. Scr. 76, 146 (2007)
%%Wu and Yang, PRA 56, 2443 (1997) -- though it should be noted that
%%the diagonalization used in this function is different from the paper as
%%it is in terms of qubit detunings instead of level energies.
%%  [CHI0, CHI1] = NUMERICALCAVITYMLS(F_C, F_10, ALPHA, CHI_LP, EPSILON, F_D, KAPPA, NMAX, PLOTS)
%%  Calculate the chi shifts for a two-level qubit-cavity system. All frequencies
%%  should be entered in Hz. F_C is the cavity frequency, F_10 the qubit transition
%%  frequency, ALPHA is the anharmonicity, CHI_LP the low power chi shift (used to calculate g_0), EPSILON the
%%  drive strength, f_d the strength of the cavity drive, and kappa the cavity photon
%%  number decay rate. NMAX controls the maximum cavity occupation that will be considered
%% (1e5 is a good default), while PLOTS should be a boolean that controls display of plots.
%%  If epsilon is a zero-length vector, cavity occupation as a function of drive will not
%%  be calculated.

w_c = 2*pi*f_c/1e9;
delta = 2*pi*(f_10-f_c)/1e9;
g0 = sqrt(abs(delta*1e9)*chi_lp*2*pi)/1e9;
d10 = 2*pi*(f_10-f_c)/1e9;
d21 = 2*pi*(f_10+alpha-f_c)/1e9;

E0 = zeros(Nmax,1);
E1 = zeros(Nmax,1);
E2 = zeros(Nmax,1);

for j=1:Nmax
    [E0(j), E1(j), E2(j)] = H(j-1, w_c, d10, d21, g0, sqrt(2)*g0);
end

w0 = diff(E0);
w1 = diff(E1);
w2 = diff(E2);

chi0 = (1e9*w0/2/pi - f_c);
chi1 = (1e9*w1/2/pi - f_c);
chi2 = (1e9*w2/2/pi - f_c);

if plots
    figure(1);clf;
    semilogx(1:Nmax-1, chi0/1e6, 'b-', 'LineWidth', 4); hold on;
    semilogx(2:Nmax-1, chi1(2:end)/1e6, 'r-', 'LineWidth', 4);
    semilogx(3:Nmax-1, chi2(3:end)/1e6, 'g-', 'LineWidth', 4); hold off;
    set(gca, 'FontSize', 16); grid on;
    xlabel('n', 'FontSize', 18);
    ylabel('\chi (MHz)', 'FontSize', 15);
    legend('\chi_0', '\chi_1', '\chi_2');
end

w_d = 2*pi*(f_d);

N = (1:Nmax-1)';

if ~isempty(epsilon)
    
        eps = epsilon;

    for j=1:length(eps)

        P = 4*eps(j)^2./( 4*(w_d - (2*pi*f_c - 2*pi*chi1)).^2 + (kappa)^2 );
        M = 4*eps(j)^2./( 4*(w_d - (2*pi*f_c - 2*pi*chi0)).^2 + (kappa)^2 );

        [~,ind] = min(abs(N - P));
        np(j) = N(ind);
        [~,ind] = min(abs(N - M));
        nm(j) = N(ind);

    end

    figure(3);clf;
    loglog(eps/1e6, np, 'r-', 'LineWidth', 4); hold on;
    loglog(eps/1e6, nm, 'b-', 'LineWidth', 4);
    set(gca, 'FontSize', 16); grid on;
    xlabel('\epsilon (MHz)', 'FontSize', 18);
    ylabel('n_{cav}', 'FontSize', 18)
    legend('n_1', 'n_0')
    xlim([min(eps), max(eps)]/1e6);


end
end
function [E0, E1, E2] = H(n, w_c, d10, d21, g0, g1)
%diagonalization of the hamiltonian in each block

V = [[0, g0*sqrt(n), 0]; [g0*sqrt(n), d10-w_c, g1*sqrt(n-1)]; ...
     [0, g1*sqrt(n-1), d10+d21-w_c]];

lambda = eig(V);
 
if n==0 
    E0 = w_c + lambda(1);
    E1 = 0;
    E2 = 0;
elseif n==1
    E0 = 2*w_c + lambda(1);
    E1 = w_c + lambda(2);
    E2 = 0;
else
    E0 = n*w_c + (w_c + lambda(1));
    E1 = (n-1)*w_c + w_c + lambda(2);
    E2 = (n-2)*w_c + w_c + lambda(3);
end
end

    