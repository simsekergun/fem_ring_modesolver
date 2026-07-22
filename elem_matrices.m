% =========================================================================
%  ELEMENT MATRICES  —  implements the verified 6x6 formulas exactly
% =========================================================================
function [Ke, Me] = elem_matrices(re, ze, eps_r, m, sgns, gw, gL, lp)
%
%  re, ze   : node coordinates [3x1]
%  eps_r    : relative permittivity (scalar)
%  m        : azimuthal mode number
%  sgns     : edge orientation signs [1x3], each +1 or -1
%  gw, gL   : quadrature weights [nqx1] and barycentric coords [nqx3]
%  lp       : local edge node pairs [3x2] = [1 2; 2 3; 3 1]
%
%  Returns 6x6 Ke (complex Hermitian) and Me (real symmetric).
%  DOF order: [u1 u2 u3  p1 p2 p3]  (edge, then nodal)
A = abs(polyarea(re, ze));

% ── Barycentric derivatives ──────────────────────────────────────────────
%   L_i = [const + rho*(z_j - z_k) + z*(rho_k - rho_j)] / (2A)
%   =>  dL_i/drho = (z_j - z_k)/(2A),   dL_i/dz = (rho_k - rho_j)/(2A)
c = [ze(2)-ze(3); ze(3)-ze(1); ze(1)-ze(2)] / (2*A);   % dL/drho
b = [re(3)-re(2); re(1)-re(3); re(2)-re(1)] / (2*A);   % dL/dz

% ── 2D scalar curl of each Nedelec basis (constant per element) ──────────
%   gamma_k = s_k * 2*(b_ni*c_nj - b_nj*c_ni)
gamma = zeros(3,1);
for k = 1:3
    ni = lp(k,1);  nj = lp(k,2);
    gamma(k) = sgns(k) * 2*(b(ni)*c(nj) - b(nj)*c(ni));
end
Ke = zeros(6,6);   % complex
Me = zeros(6,6);   % real
nq = length(gw);
for q = 1:nq
    L    = gL(q,:)';               % barycentric coords [3x1]
    rhoq = L' * re;                % physical rho at quad point
    % Quadrature weights:
    %   wt_r = gw*2A*rho  ->  for K_ee, K_pp, M_ee, M_pp
    %   wt_1 = gw*2A      ->  for K_ep, K_pe  (rho cancels 1/rho)
    wt_r = gw(q) * 2*A * rhoq;
    wt_1 = gw(q) * 2*A;

    % ── Nedelec basis vectors at quad point ─────────────────────────────
    %   N_k = s_k*(L_ni*gradL_nj - L_nj*gradL_ni)
    %   gradL_i = [c_i; b_i]   (rho-component first, z-component second)
    Nrho = zeros(3,1);   % rho-component of each Nedelec basis
    Nz   = zeros(3,1);   % z-component
    for k = 1:3
        ni = lp(k,1);  
        nj = lp(k,2);
        Nrho(k) = sgns(k) * (L(ni)*c(nj) - L(nj)*c(ni));
        Nz(k)   = sgns(k) * (L(ni)*b(nj) - L(nj)*b(ni));
    end
    % ── Loop over DOF pairs ──────────────────────────────────────────────
    for i = 1:3
        for j = 1:3
            % ── K_ee(i,j) ────────────────────────────────────────────────
            %   = INT{ gamma_i*gamma_j + (m/rho)^2*(N_i.N_j) } * rho dA
            Ke(i,j) = Ke(i,j) + wt_r * ( gamma(i)*gamma(j) ...
                        + (m/rhoq)^2 * (Nrho(i)*Nrho(j) + Nz(i)*Nz(j)) );
            % ── K_ep(i,j) ────────────────────────────────────────────────
            %   = INT{ jm*[ b_j*N_i^z + (c_j + L_j/rho)*N_i^rho ] } * dA
            Ke(i, 3+j) = Ke(i, 3+j) + wt_1 * 1j*m * ...
                          ( b(j)*Nz(i) + (c(j) + L(j)/rhoq)*Nrho(i) );
            % ── K_pp(i,j) ────────────────────────────────────────────────
            %   = INT{ b_i*b_j + (c_i+L_i/rho)*(c_j+L_j/rho) } * rho dA
            ci_rho = c(i) + L(i)/rhoq;
            cj_rho = c(j) + L(j)/rhoq;
            Ke(3+i, 3+j) = Ke(3+i, 3+j) + wt_r * ( b(i)*b(j) + ci_rho*cj_rho );
            % ── M_ee(i,j) ────────────────────────────────────────────────
            Me(i,j) = Me(i,j) + wt_r * eps_r * (Nrho(i)*Nrho(j) + Nz(i)*Nz(j));
            % ── M_pp(i,j) ────────────────────────────────────────────────
            Me(3+i, 3+j) = Me(3+i, 3+j) + wt_r * eps_r * L(i)*L(j);
        end
    end
end
% K_pe = K_ep^H 
Ke(4:6, 1:3) = Ke(1:3, 4:6)';
end