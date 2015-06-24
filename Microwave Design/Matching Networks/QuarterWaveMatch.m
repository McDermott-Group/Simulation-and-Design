

%Matching to a QPC with a high-impedance quarter wave line


%frequency range
f = 1e9:1e7:20e9;

%wire bonds
L1 = 0e-9;          %board to chip
L2 = 1e-9;          %chip to qpc

Cg = 150e-15;       %ground cap to qpc

Rqpc = 20e3;        %qpc impedance

ZL = 1000;         %transmission line 
f0 = 5.3e9;

Z0 = 50;

%%%%%%%

om = 2*pi*f;
bl = pi*f/2/f0;


Z1 = 1i*om*L2 + 1./(1i*om*Cg + 1/Rqpc);

Z1 = 1./(1./Z1 + 1./(1i*om*8e-9));

Ztl = ZL*(Z1 + 1i*ZL*tan(bl))./(ZL + 1i*Z1.*tan(bl));

Zin = 1i*om*L1 + Ztl;


Gamma = (Zin-Z0)./(Zin+Z0);

SWR = (1+abs(Gamma))./(1-abs(Gamma));

[FF, CC] = meshgrid(f, ZL);

figure(1);
semilogy(f/1e9, SWR, 'k-', 'LineWidth', 2);
grid on;
set(gca, 'FontSize', 14);
xlabel('Frequency [GHz]', 'FontSize', 14);
ylabel('SWR, Z_0 = 50 \Omega', 'FontSize', 14);

figure(2);
plot(f/1e9, real(Z1), 'r-', 'LineWidth', 2);
grid on;
set(gca, 'FontSize', 14);
xlabel('Frequency [GHz]', 'FontSize', 14);
ylabel('Re Z_{qpc}', 'FontSize', 14);


figure(3);
plot(f/1e9, imag(Z1), 'r-', 'LineWidth', 2);
grid on;
set(gca, 'FontSize', 14);
xlabel('Frequency [GHz]', 'FontSize', 14);
ylabel('Im Z_{qpc}', 'FontSize', 14);

