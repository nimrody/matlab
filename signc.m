function y = signc(x)

% y = signc(x)
%
%   y = sign(real(x)) + 1j*sign(imag(x))

y = sign(real(x)) + 1j*sign(imag(x));