function [ap, am] = cavitySimulationChiShift(time, drive, f_dr, n, chi0, chi1, f_cav, kappa, plots)
%CAVITYSIMULATIONCHISHIFT Simulate the trajectory of a linear cavity with
%an amplitude-dependent chi shift.
%see for example Gambetta et al. PRA 74, 042318 (2006)
ap0 = 0;
am0 = 0;

aplus = @(t,a) cavity(a, drive, f_dr, f_cav, n, chi1, kappa);
amin = @(t,a) cavity(a, drive, f_dr, f_cav, n, chi0, kappa);
[~,ap] = ode45(aplus, time, ap0);
[~,am] = ode45(amin, time, am0);
    
if plots
    figure(12);clf;
    hold on;
    plot(imag(ap), real(ap), 'r-', 'LineWidth', 2);
    plot(imag(am), real(am), 'b-', 'LineWidth', 2)
    xlabel('Im |\alpha_\pm|', 'FontSize', 15);
    ylabel('Re |\alpha_\pm|', 'FontSize', 15);
    set(gca, 'FontSize', 14); grid on;
    legend('\alpha_+', '\alpha_-')
    axis equal;
    
    figure(23);clf;
    hold on;
    plot3(imag(ap), real(ap), time/max(time), 'r-', 'LineWidth', 2);
    plot3(imag(am), real(am), time/max(time),'b-', 'LineWidth', 2)
    xlabel('Im |\alpha_\pm|', 'FontSize', 15);
    ylabel('Re |\alpha_\pm|', 'FontSize', 15);
    set(gca, 'FontSize', 14); grid on;
    legend('\alpha_+', '\alpha_-')
    
    figure(28); clf;
    hold on;
    plot(time, abs(ap).^2, 'r-', 'LineWidth', 2);
    plot(time, abs(am).^2, 'b-', 'LineWidth', 2);
    xlabel('Time (ns)', 'FontSize', 15);
    ylabel('|\alpha_\pm|^2', 'FontSize', 15);
    set(gca, 'FontSize', 14); grid on;
    legend('\alpha_+', '\alpha_-')
    
    n_d = (2*drive/kappa)^2
    
    figure(49); clf;
    plot(time, abs(ap-am).^2/n_d, 'k-', 'LineWidth', 2);
    xlabel('Time (ns)', 'FontSize', 15);
    ylabel('|\alpha_+ - \alpha_-|^2/n_{dr} (MHz)', 'FontSize', 15);
    set(gca, 'FontSize', 14); grid on;

end
end

function adot = cavity(alpha, epsilon, f_dr, f_cav, n, chi, kappa)

    Ncav = abs(alpha).^2;
    if Ncav <= 2
        Chi = chi(1);
    elseif Ncav > n(end)
        Chi = 0;
    else
        Chi = interp1(n, chi, Ncav);
    end
    
    adot = 1i*epsilon/1e9 - 1i*(2*pi*(f_dr - (f_cav - Chi))/1e9 - 1i*kappa/2e9)*alpha;
    
end

