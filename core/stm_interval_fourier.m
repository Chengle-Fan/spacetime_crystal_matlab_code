function value = stm_interval_fourier(n, xa, xb, Lambda)
%STM_INTERVAL_FOURIER Fourier integral over one spatial subinterval.
%
% value = (1/Lambda)*integral_xa^xb exp(-i*n*g*x) dx, g=2*pi/Lambda.

if n == 0
    value = (xb - xa)/Lambda;
else
    g = 2*pi/Lambda;
    value = (exp(-1i*n*g*xb) - exp(-1i*n*g*xa)) ...
        /(-1i*n*g*Lambda);
end
end
