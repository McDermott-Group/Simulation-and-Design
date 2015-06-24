%CoupledCPW.m
%Guilhem Ribeill, 6/8/2011
%McDermott Group, UW-Madison
%
%This file calculates the physical parameters of edge-coupled CPW lines.
%Equations are based off of "Coplanar Waveguide Circuits, Components, and
%Systems", Rainee Simons, chapter 7.6

%INPUTS:
%e_r : relative permittivity of substrate
%h   : thickness of substrate, microns
%S   : strip width, microns
%C   : voltage coupling coefficient (dB)
%Z0  : feedline impedance
%f0  ; Center frequency, Hz
%Ld  : Length of the coupled lines in degrees

function CoupledCPW(e_r, h, S, CdB, Z0, f0)

    %Maxmimum number of iterations
    MAXITER = 100000;
    %Tolerance for gradient-descent algorithm
    tol = 1e-5;
    %Gradient descent step size
    gamma = 0.01;

    %calcualte odd and even mode impedances
    C = 10^(-CdB/20);
    
    Z0e_r = Z0*(4*C - 3 + sqrt(9 - 8*C^2))/(2*C*sqrt((1-C)/(1+C)));
    Z0o_r = Z0*(4*C + 3 - sqrt(9 - 8*C^2))/(2*C*sqrt((1+C)/(1-C)));
    
    %initial guesses
    d0 = 0.5*S;
    W0 = 1.25*S;
    
    %solution vector
    x = [d0; W0];
    
    G = zeros(2,1);
    J = zeros(2,2);
    
    for step=1:MAXITER
        
        G(1,1) = EvenMode(x(1), S, x(2), e_r, h) - Z0e_r;
        G(2,1) = OddMode(x(1), S, x(2), e_r, h)  - Z0o_r;
        
        t(step) = abs(0.5*G'*G);
        
        if (t(step) < tol)
            break;
        end
        
        J = Jacobian(x, S, e_r, h);
        
        x = x - gamma*J'*G;
    end
    
    %plot the target function vs. step
    semilogy(1:step, t, 'k.');
    grid on;
    set(gca, 'FontSize', 16);
    xlabel('Step #', 'FontSize', 16);
    ylabel('F(x)', 'FontSize', 16);

    
    %if the simulation has failed, exit
    if (step == MAXITER)
        disp('FAILED TO CONVERGE...');
        return;
    end
    
    
    d = 6;%x(1);
    W = 74.5;%x(2);
    
    %Output Relevant Data
    [Z0e, e_e] = EvenMode(d, S, W, e_r, h);
    [Z0o, e_o] = OddMode(d, S, W, e_r, h);
    
    [Lopt, Ceff] = OptimalLength(Z0e, Z0o, e_e, e_o, Z0, f0);
    
    disp(' ');
    disp('Coupled CPW Calculator Results:');
    disp(' ');
    
    disp(['Target coupling: ' num2str(CdB) ' dB']);
    disp(['Calc. coupling: ' num2str(-Ceff) ' dB']);
    
    disp(['Target Z0e: ' num2str(Z0e_r) ' Ohms']);
    disp(['Target Z0e: ' num2str(Z0e) ' Ohms']);
    disp(['Calc Z0o: ' num2str(Z0o_r) ' Ohms']);
    disp(['Calc Z0o: ' num2str(Z0o) ' Ohms ']);
    
    disp(['Even Mode Phase Velcity: ', num2str(1/sqrt(e_e)), ' c']);
    disp(['Odd Mode Phase Velcity: ', num2str(1/sqrt(e_o)), ' c']);
    
    disp(['Conductor Width: ' num2str(S) ' microns']);
    disp(['Gap Width: ' num2str(W) ' microns']);
    disp(['Conductor Separation: ' num2str(d) 'microns']);
    disp(['Total Length: ' num2str(Lopt/1e-6) ' microns']);
 

end

%calculate the coupling coefficient vs. line length
function [Lopt, Ceff] = OptimalLength(Z0e, Z0o, e_e, e_o, Z0, f0)

    %speed of light
    c = 3e8;
    %lambda/4 length (maximum coupler length)
    Lmax = c/4/f0;
    dL = 1e-8; %step length
    N = Lmax/dL;
    L = linspace(0, Lmax, N);
    
    %even and odd electrical lengths
    theta_e = 2*pi*L*f0*sqrt(e_e)/c;
    theta_o = 2*pi*L*f0*sqrt(e_o)/c;
    
    %input impedances
    Zi_e = Z0e*( sqrt(Z0o) + 1i*tan(theta_e)*sqrt(Z0e) )./( sqrt(Z0e) + 1i*tan(theta_e)*sqrt(Z0o) );
    Zi_o = Z0o*( sqrt(Z0e) + 1i*tan(theta_o)*sqrt(Z0o) )./( sqrt(Z0o) + 1i*tan(theta_o)*sqrt(Z0e) );
    
    %coupling coeff
    C = Zi_e./(Zi_e+Z0) - Zi_o./(Zi_o+Z0);
    
    figure(2);
    plot(L/1e-6, abs(C).^2, 'r.');
    grid on;
    set(gca, 'FontSize', 16);
    xlabel('Coupler Length (\mu m)', 'FontSize', 16);
    ylabel('Coupling (dB)', 'FontSize', 16);
    
    [Ceff, Cind] = max(abs(C).^2);
    
    Ceff = 10*log10(Ceff);
    Lopt = L(Cind);
end
    
    
    

%calculate the jacobian matrix
function J = Jacobian(x, S, e_r, h)

    W = x(2);
    d = x(1);

    dh = sqrt(eps)*d;
    temp = d+dh;
    dh = temp-d;
    
    Z0e = EvenMode(d, S, W, e_r, h);
    Z0o = OddMode(d, S, W, e_r, h);
    
    J(1,1) = ( EvenMode(d+dh, S, W, e_r, h) - Z0e )/dh;
    J(2,1) = ( OddMode(d+dh, S, W, e_r, h) - Z0o )/dh;
    
    dh = sqrt(eps)*W;
    temp = W+dh;
    dh = temp-W;
    
    J(1,2) = ( EvenMode(d, S, W+dh, e_r, h) - Z0e )/dh;
    J(2,2) = ( OddMode(d, S, W+dh, e_r, h) - Z0o )/dh;
    
end
    

function [Zoe, e_e] = EvenMode(d, S, W, e_r, h)

    r  = d/(d+S);
    k1 = (d+2*S)/(d+2*S+2*W);
    delta = sqrt( (1-r^2)/(1-k1^2*r^2) );
    
    r1 = sinh(pi*d/4/h)/sinh( (pi/2/h)*(d/2+S) );
    k2 = sinh( (pi/2/h)*(d/2+S) )/sinh( (pi/2/h)*(d/2+S+W) );
    
    psi = sqrt( (1 - r1^2)/(1 - k2^2*r1^2) );
    
    e_e = 1 + 0.5*(e_r-1)*K(psi*k2)*Kp(delta*k1)/Kp(psi*k2)/K(delta*k1);
    
    Zoe = 60*pi*Kp(delta*k1)/sqrt(e_e)/K(delta*k1);
    
end

function [Zoo, e_o] = OddMode(d, S, W, e_r, h)

    r  = d/(d+S);
    k1 = (d+2*S)/(d+2*S+2*W);
    delta = sqrt( (1-r^2)/(1-k1^2*r^2) );

    C13 = sinh( (pi/2/h)*(d/2+S) )^2;
    C14 = sinh(pi*d/4/h)^2;
    C15 = sinh((pi/2/h)*(d/2+S+W))^2;
    
    C11 = 0.5*( ((1+C13)/(1+C14))^(0.25) - ((1+C13)/(1+C14))^(-0.25) );
    C12 = 0.5*( sqrt(1 + C15)/(((1+C13)*(1+C14))^(0.25)) - (((1+C13)*(1+C14))^(0.25))/sqrt(1 + C15));
    
    chi = -0.5*( ((1+C13)*(1+C14))^(0.25) - ((1+C13)*(1+C14))^(-0.25) );
    
    kappa = (1/(C12 - chi))*(-1 - C12*chi/C11^2 - ( (C12^2/C11^2 - 1)*(chi^2/C11^2-1))^0.5);
    
    k3 = C11*(1 + kappa*C12)/(C12 + kappa*C11^2);
    
    e_o = 1 + (e_r-1)*K(k3)*Kp(delta)/Kp(k3)/K(delta);
    
    Zoo = 60*pi*Kp(delta)/sqrt(e_o)/K(delta);
    
end
    
    

%complete elliptic integral of the first kind
function out = K(x)

    out = ellipke(x^2);
    
end

%derivative of the complete elliptic integral of the first kind
function out = Kp(x)

    out = ellipke(1-x^2);
end

