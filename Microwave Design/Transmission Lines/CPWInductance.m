%Calculate inductance and other parameters of CPW line including
%correction for kinetic inductance

function CPWInductance(e_r, h, S, W, t, lambda, f0)

    e0 = 8.854e-12;
    m0 = 4*pi*1e-7;

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
    
    L = 1./(4*f0*1e9*sqrt(Ltot*C));
    
    disp(['Strip Width: ' num2str(S) ' microns']);
    disp(['Gap Width:   ' num2str(W) ' microns']);
    disp(' ');
    
    disp(['Impedance: ' num2str(Z) ' ohms']);
    disp(['Capacitance per unit length: ' num2str(C*1e6) ' pF/micron']);
    disp(['Geometric inductance per unit length: ' num2str(Lg*1e6) ' pH/micron']);
    disp(['Kinetic inductance per unit length: ' num2str(Lk*1e6) ' pH/micron']);
    disp(['Total inductance per unit length: ' num2str(Ltot*1e6) ' pH/micron']);
    disp(['Length of a quarter-wave section at ' num2str(f0) ' GHz: ' [num2str(L/1e-6)] ' microns']);
    
end

%complete elliptic integral of the first kind
function out = K(x)

    out = ellipke(x^2);
    
end

%derivative of the complete elliptic integral of the first kind
function out = Kp(x)

    out = ellipke(1-x^2);
end