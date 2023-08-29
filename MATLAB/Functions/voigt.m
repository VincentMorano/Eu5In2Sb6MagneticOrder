function voiProfile = voigt(x, x0, sig, gam)
%voigt Voigt function
%   Given the desired x positions (e.g. a3) and parameters for the center
%   x0, sigma, and gamma, return the voigt function values at each x. A
%   prefactor is not included, this is simply the convolution of a Gaussian
%   and a Lorentzian function.

    voiProfile = zeros(length(x), 1);
    gau = @(var, sig) exp(-var.^2./(2.*sig.^2))./(sig.*sqrt(2*pi));
    lor = @(var, gam) gam./pi./(var.^2+gam.^2);
    for i = 1:length(x)
        fun = @(var) gau(var, sig).*lor((x(i) - x0)-var, gam);
        voiProfile(i) = integral(fun, -inf, inf);
    end
end