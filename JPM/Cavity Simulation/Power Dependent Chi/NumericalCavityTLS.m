function [chi0, chi1] = NumericalCavityTLS(f_c, f_10, chi_lp, epsilon, ...
                                            f_d, kappa, Nmax, plots)
%%NUMERICALCAVITYTLS -- Calculate chi shifts and cavity state as a function
%%of drive for a two level qubit. The cavity is treated semiclasically, and
%%the qubit is treated using an exact diagonalization of the
%%Jaynes-Cummings hamiltonian in the total excitation number basis. See:
%%Boissoneault et al., PRL 105, 100504 (2010)
%%Bishop et al., PRL 105, 100505 (2010)
%%Larson and Stenholm, PRA 73, 033805 (2006)
%%  [CHI0, CHI1] = NUMERICALCAVITYTLS(F_C, F_10, CHI_LP, EPSILON, F_D, KAPPA, NMAX, PLOTS)
%%  Calculate the chi shifts for a two-level qubit-cavity system. All frequencies
%%  should be entered in Hz. F_C is the cavity frequency, F_10 the qubit transition
%%  frequency, CHI_LP the low power chi shift (used to calculate g_0), EPSILON the
%%  drive strength, f_d the strength of the cavity drive, and kappa the cavity photon
%%  number decay rate. NMAX controls the maximum cavity occupation that will be considered
%% (1e5 is a good default), while PLOTS should be a boolean that controls display of plots.
%%  If epsilon is a zero-length vector, cavity occupation as a function of drive will not
%%  be calculated.


w_c = 2 * pi * f_c / 1e9;
delta = 2 * pi * (f_10 - f_c) / 1e9;
g = sqrt(abs(delta * 1e9) * chi_lp * 2 * pi) / 1e9;

Ep = zeros(Nmax,1);
Em = zeros(Nmax,1);

for j=1:Nmax
    [Ep(j), Em(j)] = H(j - 1, w_c, delta, g);
end

wp = diff(Ep);
wm = diff(Em);

chi0 = (1e9*wp/2/pi - f_c);
chi1 = (1e9*wm/2/pi - f_c);

if plots

    figure(1);clf;
    semilogx(0:Nmax-2, (wp/2/pi - f_c/1e9)*1e3, 'r-', 'LineWidth', 4); hold on;
    semilogx(0:Nmax-2, (wm/2/pi - f_c/1e9)*1e3, 'b-', 'LineWidth', 4); hold off;
    set(gca, 'FontSize', 16); grid on;
    xlabel('n', 'FontSize', 18);
    ylabel('\chi (MHz)', 'FontSize', 18)
    legend('\chi_1', '\chi_0')
    
end

if ~isempty(epsilon) 
    eps = epsilon;
    w_d = 2*pi*(f_d);

    N = (1:Nmax-1)';
    
    np = zeros(length(eps),1);
    nm = zeros(length(eps),1);
    
    for j=1:length(eps)

        P = 4 * eps(j)^2 ./ ...
            (4 * (d - (2 * pi * f_c - 2 * pi * chi1)).^2 + kappa^2);
        M = 4 * eps(j)^2 ./ ...
            (4 * (w_d - (2 * pi * f_c - 2 * pi * chi0)).^2 + kappa^2);

        [~,ind] = min(abs(N - P));
        np(j) = N(ind);
        [~,ind] = min(abs(N - M));
        nm(j) = N(ind);

    end
    
    if plots
        figure(3);clf;
        semilogy(20*log10(eps/1e6), np, 'r-', 'LineWidth', 4); hold on;
        semilogy(20*log10(eps/1e6), nm, 'b.', 'MarkerSize', 8);
        set(gca, 'FontSize', 16); grid on;
        xlabel('\epsilon (MHz)', 'FontSize', 18);
        ylabel('n_{cav}', 'FontSize', 18)
        legend('n_1', 'n_0')
    end
end
end

function [Ep, Em] = H(n, w_c, delta, g)
%diagonalization of the hamiltonian for each block.

if n==0
    Ep = 0;
    Em = -abs(delta)/2;
else
    q = sqrt(delta^2 + 4*g^2*n);
    Ep = (n-1)*w_c + 0.5*(w_c-q);
    Em = n*w_c - 0.5*(w_c-q);
end

end