
%calculate lengths for a single stub tuner matched to 50 ohms
clc;clear all;

Z0 = 7;

f = 5.3e9;

ZL = 1i*f*2*pi*1e-9 + 1/(1i*2*pi*f*200e-15 + 1/20e3);

RL = real(ZL);
XL = imag(ZL);

t1 = (XL + sqrt(RL*((Z0-RL)^2 + XL^2)/Z0))/(RL - Z0);
t2 = (XL - sqrt(RL*((Z0-RL)^2 + XL^2)/Z0))/(RL - Z0);

if (t1 >= 0)
    
    d1 = atan(t1)/2/pi;
else
    
    d1 = (pi + atan(t1))/2/pi;    
end

if (t2 >= 0)
    
    d2 = atan(t2)/2/pi;
else
    
    d2 = (pi + atan(t2))/2/pi;    
end

B1 = (RL^2*t1 - (Z0-XL*t1)*(XL+Z0*t1))/(Z0*(RL^2 + (XL+Z0*t1)^2));
B2 = (RL^2*t2 - (Z0-XL*t2)*(XL+Z0*t2))/(Z0*(RL^2 + (XL+Z0*t2)^2));

l1 = atan(1/(B1*Z0))/2/pi;
l2 = atan(1/(B2*Z0))/2/pi;

if (l1 < 0)
    l1 = l1+0.5;
end
if (l2 < 0)
    l2 = l2+0.5;
end

disp('Single Stub Tuner Calculation');
disp('--------');
disp(['Match to Load: ' num2str(ZL) ' Ohms at ' num2str(f/1e9) ' GHz']);
disp(' ');
disp('Solution #1: ')
disp(['Distance to short-circuit stub: ' num2str(d1) ' wavelengths']);
disp(['Length of short-circuit stub: ' num2str(l1) ' wavelengths']);
disp(' ');
disp('Solution #2: ')
disp(['Distance to short-circuit stub: ' num2str(d2) ' wavelengths']);
disp(['Length of short-circuit stub: ' num2str(l2) ' wavelengths']);