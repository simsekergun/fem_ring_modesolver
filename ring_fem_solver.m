%% ── Constants ─────────────────────────────────────────────────────────────
c0  = 299792458;
mu0 = 4*pi*1e-7;
eps0 = 1/(mu0*c0^2);

%% refractive indices
n_ring     = get_refractive_index(ring_material,lambda0*1e6);
n_clad     = get_refractive_index(clad_material,lambda0*1e6);
k0         = 2*pi/lambda0; % wavenumber

rho_inner   = R_central - ring_width/2;
rho_outer  = R_central + ring_width/2;
rho_min  = rho_inner - pad_rho;
rho_max  = rho_outer + pad_rho;
z_min    = -ring_height/2 - pad_z;
z_max    =  ring_height/2 + pad_z;

d_r = lambda0/meshEPW/n_ring;
d_c = lambda0/meshEPW/n_clad;

Nr_ring = round(ring_width/d_r);
Nz_ring = round(ring_height/d_r);
Nr_clad  = round((rho_max-rho_min)/d_c);
Nz_clad = round((z_max-z_min)/d_c);

disp(['Number of elements: ' int2str([Nr_ring Nz_ring Nr_clad Nz_clad])]);

%% ── Azimuthal mode number ─────────────────────────────────────────────────
% Phase-matching: m = round(n_eff * k0 * R)
a_value = sqrt(ring_width*ring_height)/pi;
V = (2*pi*a_value/lambda0) * sqrt(n_ring^2 - n_clad^2);
n_eff_est = n_clad +(n_ring- n_clad)*(1-exp(-2*V/pi)) - lambda0/(4*pi*R_central);
m_az = round(n_eff_est * k0 * R_central);

fprintf('Azimuthal mode number  m = %d\n', m_az);

%%  Generate mesh
[rho_all, z_all, elements, eps_elem] = generate_mesh( ...
    rho_inner, rho_outer, ring_height, ...
    rho_min,   rho_max,   z_min, z_max, ...
    n_ring,    n_clad, ...
    Nr_ring,   Nz_ring, ...
    Nr_clad,   Nz_clad);

numNodes = size(rho_all, 1);
numElem  = size(elements, 1);
fprintf('Nodes=%d   Elements=%d\n', numNodes, numElem);


if mesh_plot ==1
    plot_mesh(rho_all, z_all, elements, eps_elem, rho_inner, rho_outer, ring_height, save_figures, filename)
end

filenamemode = [filename int2str(m_az)];


%% ── Edge table ────────────────────────────────────────────────────────────
[edges, elemEdges, edgeSigns] = build_edges(elements);
numEdges = size(edges,1);
fprintf('Edges=%d\n', numEdges);

Ndof_e = numEdges;
Ndof_p = numNodes;
Ndof   = Ndof_e + Ndof_p;

%% ── Quadrature ────────────────────────────────────────────────────────────
[gw, gL] = gauss7();   % 7-point, exact for degree 5

%% ── Assembly ──────────────────────────────────────────────────────────────
fprintf('Assembling...\n');

% Triplet storage (36 entries per element)
ntrip = numElem*36;
Ki = zeros(ntrip,1);  Kj = zeros(ntrip,1);  Kv = zeros(ntrip,1,'like',1+1j);
Mi = zeros(ntrip,1);  Mj = zeros(ntrip,1);  Mv = zeros(ntrip,1);
pk = 0;  pm = 0;

lp = [1 2; 2 3; 3 1];   % local edge node pairs

for e = 1:numElem
    nd   = elements(e,:);
    eids = elemEdges(e,:);    % global edge IDs for this element
    sgns = edgeSigns(e,:);    % orientation signs (+1 or -1)

    re = rho_all(nd);
    ze = z_all(nd);
    eps_r = eps_elem(e);

    % 6 global DOFs: [edge1 edge2 edge3  phi1 phi2 phi3]
    dofs = [eids,  Ndof_e + nd];

    [Ke, Me] = elem_matrices(re, ze, eps_r, m_az, sgns, gw, gL, lp);

    for ii = 1:6
        for jj = 1:6
            pk = pk+1;  Ki(pk) = dofs(ii);  Kj(pk) = dofs(jj);  Kv(pk) = Ke(ii,jj);
            pm = pm+1;  Mi(pm) = dofs(ii);  Mj(pm) = dofs(jj);  Mv(pm) = Me(ii,jj);
        end
    end
end

K = sparse(Ki(1:pk), Kj(1:pk), Kv(1:pk), Ndof, Ndof);
M = sparse(Mi(1:pm), Mj(1:pm), Mv(1:pm), Ndof, Ndof);
K = (K + K')/2;   % enforce Hermitian
M = (M + M')/2;
fprintf('Assembly done.  nnz(K)=%d\n', nnz(K));

%% ── Boundary conditions (PEC on outer box) ────────────────────────────────
tol_bc = 1e-10 * max(rho_max-rho_min, z_max-z_min);
bNodes = find( abs(rho_all-rho_min)<tol_bc | abs(rho_all-rho_max)<tol_bc | ...
               abs(z_all  -z_min  )<tol_bc | abs(z_all  -z_max  )<tol_bc );

% Boundary edge DOFs: edges whose both nodes are on the boundary
isBN = false(numNodes,1);  isBN(bNodes) = true;
bEdges = find(isBN(edges(:,1)) & isBN(edges(:,2)));

% Boundary nodal DOFs
bPhi = Ndof_e + bNodes;

fixedDOF = unique([bEdges(:); bPhi(:)]);
% freeDOF  = setdiff((1:Ndof)', fixedDOF);

% PEC boundary condition is removed
%Kf = K(freeDOF, freeDOF);
%Mf = M(freeDOF, freeDOF);
%fprintf('Free DOFs: %d\n', length(freeDOF));
freeDOF = (1:Ndof)';
Kf = K;
Mf = M;

%% ── Shift-invert eigensolver ──────────────────────────────────────────────
% Eigenvalue = k0^2.  Shift to free-space k0^2.
sigma  = k0^2;
nev    = number_of_modes;
fprintf('Solving (sigma = %.6e m^{-2})...\n', sigma);
opts.tol   = 1e-10;
opts.maxit = 3000;
opts.isreal = false;
[V, D] = eigs(Kf, Mf, nev, sigma, opts);
lam = diag(D);
fprintf('Done.\n');

%% ── Filter and sort ───────────────────────────────────────────────────────
% Keep modes with real positive eigenvalue and small imaginary part
ok   = real(lam) > 0  &  abs(imag(lam)) < 0.05*abs(real(lam));
lam  = lam(ok);
V    = V(:, ok);
[~, si] = sort(abs(real(lam) - k0^2));   % sort by proximity to target
lam  = lam(si);
V    = V(:, si);

k0eig = sqrt(real(lam));
neff = m_az ./ (k0eig * R_central);

freqs = k0eig * c0 / (2*pi);


nplot = min(number_of_modes, numel(lam));

%% ── Reconstruct and plot modes ────────────────────────────────────────────
modes = zeros(Ndof, nplot);
for k = 1:nplot
    modes(freeDOF, k) = V(:, k);
end

dominant_mode = zeros(number_of_modes,1);

    for k = 1:nplot
        [Er, Ep, Ez, Gamma] = reconstruct(modes(:,k), rho_all, z_all, elements, ...
            eps_elem, Ndof_e, elemEdges, edgeSigns, lp,n_ring,n_clad);        
        [~, dominant_mode(k)] = max([mean(abs(Er)) mean(abs(Ep)) mean(abs(Ez))]);

        [Hr, Hp, Hz] = compute_Hfield(modes(:,k), rho_all, z_all, elements, ...
        eps_elem, Ndof_e, elemEdges, edgeSigns, lp, n_ring, n_clad, m_az, lambda0);

        %% ── Power normalization ────────────────────────────────────────────────
        % Compute unnormalized guided power P0 = (1/2) Re{ integral (Er*Hp* - Ep*Hr*) dA }
        P0 = 0;
        for e = 1:numElem
            nd  = elements(e,:);
            re  = rho_all(nd)*1e6;  % convert to um 
            ze = z_all(nd)*1e6;      % convert to um 
            A   = abs(polyarea(re, ze));
            rhoq = mean(re);                      % centroid rho (matches reconstruct's convention)

            Er_c = mean(Er(nd));   Ez_c = mean(Ez(nd));   % node-averaged fields at this element
            Hr_c = mean(Hr(nd));   Hz_c = mean(Hz(nd));

            S_phi = Er_c*conj(Hz_c) - Ez_c*conj(Hr_c);    % azimuthal Poynting component
            P0 = P0 + 0.5*real(S_phi) * A;
        end

        % Scale factor to bring guided power to 1 W (or any chosen P_target)
        P_target = 1;                 % [W]
        alpha = sqrt(P_target / P0);  % real scale factor

        Er = alpha*Er;  Ep = alpha*Ep;  Ez = alpha*Ez;
        Hr = alpha*Hr;  Hp = alpha*Hp;  Hz = alpha*Hz;        
        
        if mode_plot == 1

            plot_Emode(Er, Ep, Ez, Gamma, rho_all*1e6, z_all*1e6, elements, ...
            rho_inner*1e6, rho_outer*1e6, ring_height*1e6, k, freqs(k), neff(k),...
            save_figures, filenamemode,k);

            plot_Hmode(Hr, Hp, Hz, Gamma, rho_all*1e6, z_all*1e6, elements, ...
            rho_inner*1e6, rho_outer*1e6, ring_height*1e6, k, freqs(k), neff(k),...
            save_figures, [filenamemode 'H'],k);            
            
        end
    end

lam_res = c0./freqs*1e6;

fprintf('\n%4s  %14s %14s %10s %10s\n', 'Mode', 'lambda (um)' ,'freq [THz]', 'n_eff','Dominant');
for k = 1:nplot
    fprintf('%4d  %14.6f  %14.6f %10.5f %4d\n', k, lam_res(k), freqs(k)/1e12, neff(k),dominant_mode(k));
end
