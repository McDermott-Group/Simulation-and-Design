clc;clear all;

Z0 = 50;

f = 5.3e9;

ZL = 1i*f*2*pi*1e-9 + 1/(1i*2*pi*f*200e-15 + 1/20e3);

RL = real(ZL);
XL = imag(ZL);

ZS = Z0+1i*2*pi*f*1e-9;

RS = real(ZS);
XS = imag(ZS);

ZTL = sqrt(((RS^2 + XS^2)*RL - (RL^2+XL^2)*RS)/(RS-RL));

theta = atan(ZTL*(RL-RS)/(XS*RL - XL*RS));

disp('Transmission Line Matching Section');
disp('--------');
disp(['Load: ' num2str(ZL) ' Ohms at ' num2str(f/1e9) ' GHz']);
disp(['Source: ' num2str(ZS) ' Ohms at ' num2str(f/1e9) ' GHz']);
disp(' ');
disp(['Transmission Line Impedance: ' num2str(ZTL) ' Ohms']);
disp(['Electrical Length: ' num2str(theta) ' wavelengths' ]);