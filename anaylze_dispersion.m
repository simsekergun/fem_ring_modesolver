% this function calculates the integrated dispersion

function [mu, Dint, D1, D2, D3, D4] = anaylze_dispersion(R_central, lambdas, n_eff,central_freq)
c0  = 299792458;
omega = 2*pi*c0 ./ lambdas;
beta = n_eff .* omega / c0;
mvalues = beta*R_central;
mi = ceil(min(mvalues)):floor(max(mvalues));
pf = polyfit(mvalues, n_eff,5);
n_effi = polyval(pf, mi);

freqs = c0.*mi./n_effi/(2*pi*R_central);
% Find relative mode numbers
R = abs(freqs-central_freq);          % Find difference between pump frequency and each frequency
NN = find(R==min(R));                 % Lowest difference is the index of pump frequency
mu = mi-mi(NN);

w0 = freqs(NN);
D1 = (freqs(NN+1)-freqs(NN-1))/2;
Dint = freqs-w0-D1*mu;

pf = polyfit(mu, Dint,5);
D1 = pf(end-1);
D2 = pf(end-2)*2;
D3 = pf(end-3)*6;
D4 = pf(end-4)*24;
