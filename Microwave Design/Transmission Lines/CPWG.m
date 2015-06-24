%Calulate geometry for grounded CPW line

function CPWG(e_r, h, S, Z0, f0)

    %Maxmimum number of iterations
    MAXITER = 10000;
    %Tolerance for gradient-descent algorithm
    tol = 1e-5;
    
    W = zeros(MAXITER,1);
    
    W(1) = S/2;
    W(2) = 2*S;
    
    for step=3:MAXITER
        
        F0 = CPWCalc(e_r, h, S, W(step-2)) - Z0;
        F1 = CPWCalc(e_r, h, S, W(step-1)) - Z0;
        
        if (abs(F1) < tol)
            sol = step-1;
            break;
        end
        
        W(step) = W(step-1) - F1*(W(step-1) - W(step-2))/(F1 - F0);
        
    end
    
    if (step==MAXITER)
        return;
    end
    
    [Z, e, C_cpw] = CPWCalc(e_r, h, S, W(sol));
            
    disp('CPW Calculator Results:'); disp(' ');
    
    disp(['Target Inductance: ' num2str(Z0) ' Ohms']);
    disp(['Calc. Inductance: ' num2str(Z), ' Ohms']);
    disp(['Strip Width: ' num2str(S) ' microns']);
    disp(['Gap Width: ' num2str(W(sol)) ' microns']);
    
    vph = 1/Z/C_cpw
    L4 = vph/4/f0;
    
    disp(['Length of Quarter Wavelength Line: ' num2str(L4/1e-6) ' microns']);

    
    
    
%     disp(['Capacitance per unit length: ' num2str(C_cpw*1e6) ' pF/micron']);
%     disp(['Capacitance of ' num2str(L) ' micron line: ' num2str(C_tot*1e12) ' pF']);
        
    
    
end

%calculate impedance and effective dielectric constant
function [Z, eps, C] = CPWCalc(e_r, h, S, W)

    e0 = 8.854e-12;
    m0 = 4*pi*1e-7;
    
    t = 0.1;
    lambda = 0.1;
    

    k0 = S/(S+2*W);
    k1 = sinh(pi*S/4/h)/sinh(pi*(S+2*W)/4/h);
    
    eps = 1 + 0.5*(e_r-1)*K(k1)*Kp(k0)/Kp(k1)/K(k0);
    Z   = 30*pi*Kp(k0)/sqrt(eps)/K(k0);
    
    C = 2*(e_r-1)*K(k1)/Kp(k1) + 4*K(k0)/Kp(k0);
    C = C*e0;
    
    Lg = Z^2*C;
    
    A = -t/pi + 0.5*sqrt( (2*t/pi)^2 + S^2);
    B = S^2/4/A;
    CC = B - t/pi + sqrt((t/pi)^2+W^2);
    D = 2*t/pi + CC;
    
    Lk = m0*lambda*CC/(4*A*D*K(k0))*(1.7/sinh(t/2/lambda) + 0.4/sqrt( ((B/A)^2-1)*(1-(B/D)^2) ));
    
    Ltot = Lk + Lg;
    
    Z = sqrt(Ltot/C);
    
end

%complete elliptic integral of the first kind
function out = K(x)

    out = ellipke(x^2);
    
end

%derivative of the complete elliptic integral of the first kind
function out = Kp(x)

    out = ellipke(1-x^2);
end