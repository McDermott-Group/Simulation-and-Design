%calculate impedance and effective dielectric constant of CPW
function [Z, eps, C] = CPWCalc(e_r, h, S, W)

    e0 = 8.854e-12;
    m0 = 4*pi*1e-7;
    
    t = 0.1;
    lambda = 0.3; %usually 0.1.
    

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
    %Ltot = Lg;
    
    Z = sqrt(Ltot/C);
    Ltot
    C
    
    vph = 1/sqrt(Ltot*C) % 1/sqrt(Lg*C)
    
end

%complete elliptic integral of the first kind
function out = K(x)

    out = ellipke(x^2);
    
end

%derivative of the complete elliptic integral of the first kind
function out = Kp(x)

    out = ellipke(1-x^2);
end