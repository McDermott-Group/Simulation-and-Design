%Microstrip.m

%Synthesizes a design for a superconducting microstrip line given the
%following input variables:

%   e_r:    relative permittivity of substrate
%   h:      height of substrate (microns)
%   Z0:     desired impedance
%   f0:     center frequency, in GHz

function Microstrip(e_r, h, Z0, f0)


    %Maxmimum number of iterations
    MAXITER = 2500;
    %Tolerance for gradient-descent algorithm
    tol = 1e-3;
    
%     MSCalc(11.8, 0.15, 5)
%     MSCalc(11.8, 0.15, 12)
    
    W = zeros(MAXITER,1);
    
    W(1) = 2;
    W(2) = 20;
    
    for step=3:MAXITER
        
        F0 = MSCalc(e_r, h, W(step-2)) - Z0
        F1 = MSCalc(e_r, h, W(step-1)) - Z0
        
        if (abs(F1) < tol)
            sol = step-1;
            break;
        end
        
        W(step) = W(step-1) - F1*(W(step-1) - W(step-2))/(F1 - F0);
        
    end
    
    if (step==MAXITER)
        disp('Max iterations reached!');
        return;
    end
    
    [Z, v] = MSCalc(e_r, h, W(sol));
    
    disp('Microstrips Calculator Results:'); disp(' ');
    
    disp(['Target Inductance: ' num2str(Z0) ' Ohms']);
    disp(['Calc. Inductance: ' num2str(Z), ' Ohms']);
    disp(['Strip Width: ' num2str(W(sol)) ' microns']);
    
    disp(' ');
    
    disp(['Phase velocity: ' num2str(v) ' m/s']);
    disp(['Length of quarter-wave section: ' num2str(0.25*v/(f0*1e9)/1e-6) ' microns']);     
    
    
end

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
