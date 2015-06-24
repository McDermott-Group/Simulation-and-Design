function cs = MicrostripGap(W, h, L, er)

%from http://qucs.sourceforge.net/tech/node79.html

Q1 = 0.04598*(0.03 + (W/h)^1.23)*(0.272+0.07*er);
cs = 500*h*exp(-1.86*L/h)*Q1*(1+4.19*(1-exp(-0.785*sqrt(h/W))));
