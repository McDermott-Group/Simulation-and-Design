%interdigitated capacitor from Gupta

function Co = IDC(W, S, l, eps, n)

a = W/2;
b = (W+S)/2;

e_re = (eps+1)/2;

k = tan(a*pi/4/b)^2;

C = (e_re*10^(-3)/18/pi)*(K(k)/Kp(k))*(n-1)*l;

%disp(['Capacitance is : ' num2str(C*1000) ' fF']);

Co = C*1000;

end

%complete elliptic integral of the first kind
function out = K(x)

    out = ellipke(x^2);
    
end

%derivative of the complete elliptic integral of the first kind
function out = Kp(x)

    out = ellipke(1-x^2);
end