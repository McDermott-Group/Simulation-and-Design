%calculate impedance and effective dielectric constant
function [Z, v] = MSCalc(e_r, h, W)

    %permeability/permittivity of free space
    e0 = 8.854e-12;
    m0 = 4*pi*1e-7;
    
    t = 0.1;        %thickness of metal
    lambda = 0.1;   %penetration depth
    
    %correction factor for superconducting line
    cf =  sqrt(1 + (2*lambda/h)*coth(t/lambda));
    
    %effective dielectric constant
    eps = (e_r+1)/2 + ( (e_r-1)/2 )/sqrt(1 + 12*h/W);
    
    %line impedance (only true for W/h > 1)
    Z = 120*pi/sqrt(eps)/(W/h + 1.393 + 0.667*log(W/h + 1.444));
    
    
    Z = Z*cf;
    
    %phase velocity on line
    v = 3e8/sqrt(eps)/cf;
    
end