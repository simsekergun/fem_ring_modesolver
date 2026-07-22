function [n_effs1, n_effs2, lam_res1, lam_res2] = get_ring_neff(lambda0, ...
    R_central, ring_width, ring_height, n_ring, n_clad,rho_all, z_all, elements, eps_elem,... 
    edges, elemEdges, edgeSigns,pad_rho,pad_z,mode_plot, m_az)

c0  = 299792458;
number_of_modes = 2;

rho_inner   = R_central - ring_width/2;
rho_outer   = R_central + ring_width/2;
rho_min  = rho_inner - pad_rho;
rho_max  = rho_outer + pad_rho;
z_min    = -ring_height/2 - pad_z;
z_max    =  ring_height/2 + pad_z;

k0         = 2*pi/lambda0;

numNodes = size(rho_all, 1);
numElem  = size(elements, 1);

%%%%%%%%%
numEdges = size(edges,1);

Ndof_e = numEdges;
Ndof_p = numNodes;
Ndof   = Ndof_e + Ndof_p;

%% ── Quadrature ────────────────────────────────────────────────────────────
[gw, gL] = gauss7();   % 7-point, exact for degree 5

%% ── Assembly ──────────────────────────────────────────────────────────────

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

freeDOF = (1:Ndof)';
Kf = K;
Mf = M;


%% ── Shift-invert eigensolver ──────────────────────────────────────────────
% Eigenvalue = k0^2.  Shift to free-space k0^2.
sigma  = k0^2;
nev    = number_of_modes*2;
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
    if mode_plot == 1
        plot_mode(Er, Ep, Ez, Gamma, rho_all, z_all, elements, ...
            rho_inner, rho_outer, ring_height, k, freqs(k), neff(k),...
            save_figures, filenamemode,k);
    end
end

lam_res = c0./freqs*1e6;

fprintf('\n%4s  %14s %14s %10s %10s\n', 'Mode', 'lambda (um)' ,'freq [THz]', 'n_eff','Dominant');
for k = 1:nplot
    fprintf('%4d  %14.6f  %14.6f %10.5f %4d\n', k, lam_res(k), freqs(k)/1e12, neff(k),dominant_mode(k));
end


if dominant_mode(1) == 1        % TE
    n_effs1 = neff(1);
    lam_res1 = lam_res(1);
    n_effs2 = neff(2);
    lam_res2 = lam_res(2);
else                                            % TM
    n_effs2 = neff(1);
    lam_res2 = lam_res(1);
    n_effs1 = neff(2);
    lam_res1 = lam_res(2);
end
