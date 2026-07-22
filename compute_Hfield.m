function [Hr, Hp, Hz] = compute_Hfield(mode_vec, rho_all, z_all, elements, ...
    eps_elem, Ndof_e, elemEdges, edgeSigns, lp, n_ring, n_clad, m_az, lambda0)
% COMPUTE_HFIELD  Reconstruct H-field components from FEM eigenvector
%
% Computes H_rho, H_phi, H_z at all mesh nodes by applying the curl
% relation  curl(E) = -j*omega*mu0*H  in cylindrical coordinates with
% azimuthal dependence exp(j*m*phi).
%
% Curl equations (fields proportional to exp(j*m*phi), d/dphi -> j*m):
%   (curl E)_rho = (j*m/rho)*Ez  - dEphi/dz          = -j*omega*mu0 * Hr
%   (curl E)_phi =  dEr/dz       - dEz/drho           = -j*omega*mu0 * Hp
%   (curl E)_z   = (1/rho)*d(rho*Ephi)/drho - (j*m/rho)*Er
%                = dEphi/drho + Ephi/rho  - (j*m/rho)*Er = -j*omega*mu0 * Hz
%
% INPUTS:
%   mode_vec    - Full DOF eigenvector  [Ndof_e + Nnodes, 1]
%                 First Ndof_e entries  -> Nedelec edge DOFs (E_rho, E_z)
%                 Next  Nnodes entries  -> Lagrange node DOFs (E_phi)
%   rho_all     - [Nnodes x 1] rho-coordinates of all nodes
%   z_all       - [Nnodes x 1] z-coordinates of all nodes
%   elements    - [Nelem x 3]  node connectivity  (local nodes 1,2,3)
%   eps_elem    - [Nelem x 1]  relative permittivity per element (unused here,
%                              kept for interface consistency)
%   Ndof_e      - number of edge DOFs (= number of edges)
%   elemEdges   - [Nelem x 3]  global edge indices per element
%   edgeSigns   - [Nelem x 3]  orientation signs (+1/-1) per element edge
%   lp          - [3 x 2]      local edge-node index pairs, e.g. [1 2;2 3;3 1]
%   n_ring      - ring refractive index  (unused here, for interface compat.)
%   n_clad      - cladding refractive index (unused here)
%   m_az        - azimuthal mode number  m
%   lambda0     - free-space wavelength [m]
%
% OUTPUTS:
%   Hr  - [Nnodes x 1]  H_rho at mesh nodes  (complex, [A/m] up to normalization)
%   Hp  - [Nnodes x 1]  H_phi at mesh nodes
%   Hz  - [Nnodes x 1]  H_z   at mesh nodes
%
% ALGORITHM:
%   E fields are first reconstructed element-by-element at Gauss quadrature
%   points.  Spatial derivatives are computed analytically using the linear
%   (Lagrange / Nedelec) shape-function gradients on each triangle.
%   The resulting H-field contributions are averaged to the nodes via an
%   area-weighted accumulation.
% -------------------------------------------------------------------------

%% ── Physical constants ───────────────────────────────────────────────────
mu0 = 4*pi*1e-7;
c0  = 299792458;
k0  = 2*pi / lambda0;
omega = k0 * c0;          % angular frequency

prefactor = 1 / (-1j * omega * mu0);   % curl E = -jωμ0 H  =>  H = curl(E)/(-jωμ0)

%% ── Mesh sizes ───────────────────────────────────────────────────────────
numNodes = length(rho_all);
numElem  = size(elements, 1);

%% ── Initialise node-centred accumulators ─────────────────────────────────
Hr_acc  = zeros(numNodes, 1, 'like', 1+1j);
Hp_acc  = zeros(numNodes, 1, 'like', 1+1j);
Hz_acc  = zeros(numNodes, 1, 'like', 1+1j);
weight  = zeros(numNodes, 1);          % total area weight per node

%% ── Loop over elements ───────────────────────────────────────────────────
for e = 1:numElem

    % ----- Element geometry ----------------------------------------------
    nd  = elements(e, :);          % global node indices  [n1 n2 n3]
    re  = rho_all(nd);             % rho of nodes  [3 x 1]
    ze  = z_all(nd);               % z   of nodes  [3 x 1]

    % Signed element area  (positive when nodes are CCW)
    A2  = (re(2)-re(1))*(ze(3)-ze(1)) - (re(3)-re(1))*(ze(2)-ze(1));
    A   = A2 / 2;                  % may be negative for CW ordering
    absA = abs(A);

    % Linear shape-function gradients (constant per element)
    %   L_i on nodes:  c_i = dL_i/drho,  b_i = dL_i/dz
    %   Using the standard formula with signed area A2
    % b = zeros(3,1);  c = zeros(3,1);
    % for i = 1:3
    %     j = mod(i,   3) + 1;
    %     k = mod(i+1, 3) + 1;
    %     b(i) =  (re(j) - re(k)) / A2;   % dL_i/dz    (note sign convention)
    %     c(i) = -(ze(j) - ze(k)) / A2;   % dL_i/drho
    %     % Equivalently:
    %     %   b(i) = dL_i/dz,  c(i) = dL_i/drho
    %     % matching the main solver convention:
    %     %   c_i = (rho_k - rho_j)/(2A),  b_i = (z_j - z_k)/(2A)
    % end
    % Re-derive with the explicit convention from the main solver:
    %   local indices:  1->nd(1), 2->nd(2), 3->nd(3)
    %   c_i = (rho_k - rho_j) / (2*A)   with A = |A2|/2, signed 2A = A2
    %   b_i = (z_j  - z_k ) / (2*A)
    % This matches the above because:
    %   c_1 = (re(3)-re(2))/A2 = -(re(2)-re(3))/A2
    % Let us use the explicit convention directly to avoid confusion:
    c = zeros(3,1);  b = zeros(3,1);
    c(1) = (ze(2) - ze(3)) / A2;   b(1) = (re(3) - re(2)) / A2;
    c(2) = (ze(3) - ze(1)) / A2;   b(2) = (re(1) - re(3)) / A2;
    c(3) = (ze(1) - ze(2)) / A2;   b(3) = (re(2) - re(1)) / A2;


    % ----- DOF extraction ------------------------------------------------
    eids = elemEdges(e, :);        % global edge DOFs  [3 x 1]
    sgns = edgeSigns(e, :);        % signs             [3 x 1]

    % Edge DOF values
    e_edge = mode_vec(eids);       % [3 x 1]

    % Nodal (phi) DOF values
    e_phi_nodes = mode_vec(Ndof_e + nd);   % E_phi at nodes [3 x 1]

    % ----- Centroid evaluation of E fields -------------------------------
    % Use element centroid (barycentric coords = [1/3, 1/3, 1/3])
    % for a single representative point per element.
    L_c = [1/3; 1/3; 1/3];
    rho_c = L_c' * re;        % centroid rho
    % z_c   = L_c' * ze;      % centroid z  (not explicitly needed below)

    % E_phi at centroid (Lagrange interpolation)
    Ep_c = L_c' * e_phi_nodes;    % scalar

    % E_rho, E_z at centroid via Nedelec basis
    % N_k^rho = sgn_k * (L_ni * c_nj - L_nj * c_ni)
    % N_k^z   = sgn_k * (L_ni * b_nj - L_nj * b_ni)
    Er_c = 0;
    Ez_c = 0;
    for ked = 1:3
        ni = lp(ked, 1);   nj = lp(ked, 2);
        s  = sgns(ked);
        Nr_k = s * (L_c(ni)*c(nj) - L_c(nj)*c(ni));
        Nz_k = s * (L_c(ni)*b(nj) - L_c(nj)*b(ni));
        Er_c = Er_c + e_edge(ked) * Nr_k;
        Ez_c = Ez_c + e_edge(ked) * Nz_k;
    end

    % ----- Spatial derivatives of E fields at centroid -------------------
    %
    % dEr/dz:  Er = sum_k e_k * N_k^rho(rho,z).
    %   N_k^rho = s_k*(L_ni*c_nj - L_nj*c_ni), and L_i are linear in (rho,z)
    %   => dN_k^rho/dz = s_k*(b_ni*c_nj - b_nj*c_ni)  [constant]
    dEr_dz = 0;
    dEr_drho = 0;
    for ked = 1:3
        ni = lp(ked, 1);   nj = lp(ked, 2);
        s  = sgns(ked);
        dNr_dz   = s * (b(ni)*c(nj) - b(nj)*c(ni));
        dNr_drho = s * (c(ni)*c(nj) - c(nj)*c(ni));  % = 0 always, kept for clarity
        dEr_dz   = dEr_dz   + e_edge(ked) * dNr_dz;
        dEr_drho = dEr_drho + e_edge(ked) * dNr_drho;
    end
    % Note: dNr_drho = s*(c_ni*c_nj - c_nj*c_ni) = 0, so dEr_drho = 0
    % (the rho-component of Nedelec basis has zero drho derivative on a triangle)

    % dEz/drho: N_k^z = s_k*(L_ni*b_nj - L_nj*b_ni)
    %   dN_k^z/drho = s_k*(c_ni*b_nj - c_nj*b_ni)  [constant]
    dEz_drho = 0;
    dEz_dz   = 0;
    for ked = 1:3
        ni = lp(ked, 1);   nj = lp(ked, 2);
        s  = sgns(ked);
        dNz_drho = s * (c(ni)*b(nj) - c(nj)*b(ni));
        dNz_dz   = s * (b(ni)*b(nj) - b(nj)*b(ni));  % = 0, kept for clarity
        dEz_drho = dEz_drho + e_edge(ked) * dNz_drho;
        dEz_dz   = dEz_dz   + e_edge(ked) * dNz_dz;
    end

    % dEphi/drho and dEphi/dz (Lagrange: E_phi = sum_i e_phi_i * L_i)
    %   dEphi/drho = sum_i e_phi_i * c_i
    %   dEphi/dz   = sum_i e_phi_i * b_i
    dEp_drho = e_phi_nodes' * c;   % scalar
    dEp_dz   = e_phi_nodes' * b;   % scalar

    % ----- Curl components -----------------------------------------------
    % (curl E)_rho = jm/rho * Ez  -  dEphi/dz
    curl_r = (1j*m_az / rho_c) * Ez_c  -  dEp_dz;

    % (curl E)_phi =  dEr/dz  -  dEz/drho
    curl_p = dEr_dz  -  dEz_drho;

    % (curl E)_z   = dEphi/drho  +  Ephi/rho  -  jm/rho * Er
    curl_z = dEp_drho  +  Ep_c / rho_c  -  (1j*m_az / rho_c) * Er_c;

    % ----- H field at centroid -------------------------------------------
    Hr_c = prefactor * curl_r;
    Hp_c = prefactor * curl_p;
    Hz_c = prefactor * curl_z;

    % ----- Accumulate to nodes (area weighting) --------------------------
    for i = 1:3
        Hr_acc(nd(i))  = Hr_acc(nd(i))  + absA * Hr_c;
        Hp_acc(nd(i))  = Hp_acc(nd(i))  + absA * Hp_c;
        Hz_acc(nd(i))  = Hz_acc(nd(i))  + absA * Hz_c;
        weight(nd(i))  = weight(nd(i))  + absA;
    end

end

%% ── Normalise by total weight per node ───────────────────────────────────
nz_mask = weight > 0;
Hr = zeros(numNodes, 1, 'like', 1+1j);
Hp = zeros(numNodes, 1, 'like', 1+1j);
Hz = zeros(numNodes, 1, 'like', 1+1j);

Hr(nz_mask) = Hr_acc(nz_mask) ./ weight(nz_mask);
Hp(nz_mask) = Hp_acc(nz_mask) ./ weight(nz_mask);
Hz(nz_mask) = Hz_acc(nz_mask) ./ weight(nz_mask);
end