% =========================================================================
%  FIELD RECONSTRUCTION AT NODES
% =========================================================================
function [Erho_n, Ephi_n, Ez_n, Gamma] = reconstruct(mode, rho_all, z_all, ...
    elements, eps_elem, Ndof_e, elemEdges, edgeSigns, lp, n_ring, n_bg)

numNodes = numel(rho_all);
numElem  = size(elements,1);
eF = mode(1:Ndof_e);           % edge DOFs  -> E_rho, E_z
pF = mode(Ndof_e+1:end);       % nodal DOFs -> E_phi

Erho_sum = zeros(numNodes,1);
Ez_sum   = zeros(numNodes,1);
cnt_n    = zeros(numNodes,1);
ringE = 0;  totE = 0;

for e = 1:numElem
    nd  = elements(e,:);
    re  = rho_all(nd);  ze = z_all(nd);
    A   = abs(polyarea(re, ze));
    b   = [ze(2)-ze(3); ze(3)-ze(1); ze(1)-ze(2)] / (2*A);
    c   = [re(3)-re(2); re(1)-re(3); re(2)-re(1)] / (2*A);
    rhoq = mean(re);     % centroid
    L13  = [1/3; 1/3; 1/3];

    % Evaluate E_t = sum_k u_k * N_k at centroid
    Ev = [0; 0];
    for k = 1:3
        ni = lp(k,1);  nj = lp(k,2);
        Nk = edgeSigns(e,k) * [L13(ni)*c(nj)-L13(nj)*c(ni); ...
                                L13(ni)*b(nj)-L13(nj)*b(ni)];
        Ev = Ev + Nk * eF(elemEdges(e,k));
    end

    Ep_c = mean(pF(nd));
    en   = (abs(Ev(1))^2 + abs(Ev(2))^2 + abs(Ep_c)^2) * A * rhoq;
    totE = totE + en;
    if eps_elem(e) > (n_ring^2 + n_bg^2)/2   % rough ring test by eps
        ringE = ringE + en;
    end

    for n = 1:3
        Erho_sum(nd(n)) = Erho_sum(nd(n)) + Ev(2);
        Ez_sum(nd(n))   = Ez_sum(nd(n))   + Ev(1);
        cnt_n(nd(n))    = cnt_n(nd(n)) + 1;
    end
end

g = cnt_n > 0;
Erho_n = zeros(numNodes,1);  
Ez_n = zeros(numNodes,1);
Erho_n(g) = Erho_sum(g) ./ cnt_n(g);
Ez_n(g)   = Ez_sum(g)   ./ cnt_n(g);
Ephi_n    = pF;
Gamma     = real(ringE / max(totE, 1e-30));
end