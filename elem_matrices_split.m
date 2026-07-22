% =========================================================================
%  ELEMENT MATRICES — SPLIT BY POWER OF m  (for fixed-lambda / solve-for-m)
% =========================================================================
function [Kee_A, Kee_B, Kep_C, Kpp, Mee, Mpp] = elem_matrices_split(re, ze, eps_r, sgns, gw, gL, lp)
%  Same element geometry as elem_matrices, but returns the pieces that
%  multiply m^0, m^1 (imag), and m^2 SEPARATELY, so the global matrices
%  can be assembled once and then combined for ANY m or ANY lambda without
%  re-looping over elements.
%
%   K_ee(m)   = Kee_A + m^2 * Kee_B
%   K_ep(m)   = 1j*m * Kep_C          (K_pe = K_ep^H)
%   K_pp      = Kpp                   (no m dependence)
%   M_ee, M_pp are as before.
%
A = abs(polyarea(re, ze));
b = [ze(2)-ze(3); ze(3)-ze(1); ze(1)-ze(2)] / (2*A);
c = [re(3)-re(2); re(1)-re(3); re(2)-re(1)] / (2*A);
gamma = zeros(3,1);
for k = 1:3
    ni = lp(k,1);  nj = lp(k,2);
    gamma(k) = sgns(k) * 2*(b(ni)*c(nj) - b(nj)*c(ni));
end
Kee_A = zeros(3,3);   % gamma_i*gamma_j term            (m^0)
Kee_B = zeros(3,3);   % (1/rho^2)*(N_i.N_j) term        (coeff of m^2)
Kep_C = zeros(3,3);   % [b_j*N_i^z + (c_j+L_j/rho)*N_i^rho]  (coeff of j*m)
Kpp   = zeros(3,3);   % nodal-nodal block                (m^0)
Mee   = zeros(3,3);
Mpp   = zeros(3,3);
nq = length(gw);
for q = 1:nq
    L    = gL(q,:)';
    rhoq = L' * re;
    wt_r = gw(q) * 2*A * rhoq;
    wt_1 = gw(q) * 2*A;
    Nrho = zeros(3,1);  Nz = zeros(3,1);
    for k = 1:3
        ni = lp(k,1);  nj = lp(k,2);
        Nrho(k) = sgns(k) * (L(ni)*c(nj) - L(nj)*c(ni));
        Nz(k)   = sgns(k) * (L(ni)*b(nj) - L(nj)*b(ni));
    end
    for i = 1:3
        for j = 1:3
            Kee_A(i,j) = Kee_A(i,j) + wt_r * gamma(i)*gamma(j);
            Kee_B(i,j) = Kee_B(i,j) + wt_r * (1/rhoq^2) * (Nrho(i)*Nrho(j) + Nz(i)*Nz(j));
            Kep_C(i,j) = Kep_C(i,j) + wt_1 * ( b(j)*Nz(i) + (c(j)+L(j)/rhoq)*Nrho(i) );
            ci_rho = c(i) + L(i)/rhoq;
            cj_rho = c(j) + L(j)/rhoq;
            Kpp(i,j) = Kpp(i,j) + wt_r * ( b(i)*b(j) + ci_rho*cj_rho );
            Mee(i,j) = Mee(i,j) + wt_r * eps_r * (Nrho(i)*Nrho(j) + Nz(i)*Nz(j));
            Mpp(i,j) = Mpp(i,j) + wt_r * eps_r * L(i)*L(j);
        end
    end
end
end