k0         = 2*pi/lambda0;
lambda_target = k0^2;              % fixed eigenvalue, lambda = k0^2

n_ring     = get_refractive_index(core_material,lambda0*1e6);
n_clad     = get_refractive_index(clad_material,lambda0*1e6);

rho_inner   = R_central - ring_width/2;
rho_outer   = R_central + ring_width/2;

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
m_guess = round(n_eff_est * k0 * R_central);
fprintf('Target lambda0 = %.2f um,  guess m = %.4f\n', lambda0*1e6, m_guess);

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


%% ── Edge table & DOF bookkeeping (unchanged) ──────────────────────────────
[edges, elemEdges, edgeSigns] = build_edges(elements);
numEdges = size(edges,1);
Ndof_e = numEdges;
Ndof_p = numNodes;
Ndof   = Ndof_e + Ndof_p;
[gw, gL] = gauss7();

%% ── Assemble the SIX m-independent global blocks ONCE ──────────────────────
%fprintf('Assembling m-independent blocks...\n');
ntripEE = numElem*9;  ntripEP = numElem*9;  ntripPP = numElem*9;
[iEE,jEE,vAee,vBee,vMee] = deal(zeros(ntripEE,1));
[iEP,jEP,vCep]           = deal(zeros(ntripEP,1));
[iPP,jPP,vKpp,vMpp]      = deal(zeros(ntripPP,1));
pEE = 0; pEP = 0; pPP = 0;
for e = 1:numElem
    nd   = elements(e,:);
    eids = elemEdges(e,:);
    sgns = edgeSigns(e,:);
    re = rho_all(nd);  ze = z_all(nd);  eps_r = eps_elem(e);
    [KeeA, KeeB, KepC, Kpp_e, Mee_e, Mpp_e] = ...
        elem_matrices_split(re, ze, eps_r, sgns, gw, gL, [1 2; 2 3; 3 1]);
    for ii = 1:3
        for jj = 1:3
            pEE=pEE+1; iEE(pEE)=eids(ii); jEE(pEE)=eids(jj);
            vAee(pEE)=KeeA(ii,jj); vBee(pEE)=KeeB(ii,jj); vMee(pEE)=Mee_e(ii,jj);
            pEP=pEP+1; iEP(pEP)=eids(ii); jEP(pEP)=nd(jj);
            vCep(pEP)=KepC(ii,jj);
            pPP=pPP+1; iPP(pPP)=nd(ii);  jPP(pPP)=nd(jj);
            vKpp(pPP)=Kpp_e(ii,jj); vMpp(pPP)=Mpp_e(ii,jj);
        end
    end
end
Aee = sparse(iEE,jEE,vAee,Ndof_e,Ndof_e);   Aee = (Aee+Aee.')/2;
Bee = sparse(iEE,jEE,vBee,Ndof_e,Ndof_e);   Bee = (Bee+Bee.')/2;
Mee = sparse(iEE,jEE,vMee,Ndof_e,Ndof_e);   Mee = (Mee+Mee.')/2;
Cep = sparse(iEP,jEP,vCep,Ndof_e,Ndof_p);          % NOT symmetric — keep as is
Kpp = sparse(iPP,jPP,vKpp,Ndof_p,Ndof_p);   Kpp = (Kpp+Kpp.')/2;
Mpp = sparse(iPP,jPP,vMpp,Ndof_p,Ndof_p);   Mpp = (Mpp+Mpp.')/2;
fprintf('Done.  nnz: Aee=%d Bee=%d Cep=%d Kpp=%d\n', nnz(Aee), nnz(Bee), nnz(Cep), nnz(Kpp));

%% ── Build the quadratic-eigenvalue-problem matrices  A2 m^2 + A1 m + A0 ────
% A0 = [Aee - lambda*Mee,        0        ;   0,   Kpp - lambda*Mpp]
% A1 = [0,                    1j*Cep      ;  -1j*Cep.',      0     ]
% A2 = [Bee,                     0        ;   0,             0     ]
Z_ep = sparse(Ndof_e, Ndof_p);
Z_pp = sparse(Ndof_p, Ndof_p);
A0 = [Aee - lambda_target*Mee,  Z_ep;
    Z_ep.',                  Kpp - lambda_target*Mpp];
A1 = [sparse(Ndof_e,Ndof_e),    1j*Cep;
    -1j*Cep.',                sparse(Ndof_p,Ndof_p)];
A2 = [Bee,                      Z_ep;
    Z_ep.',                   Z_pp];


%% ── Linearize (Tisseur & Meerbergen companion form) and solve for m ───────
%   Q(m) x = (m^2 A2 + m A1 + A0) x = 0   <=>   Amat z = m Bmat z,
%   z = [m*x ; x]
I_N = speye(Ndof);
Amat = [-A1, -A0;
    I_N, sparse(Ndof,Ndof)];
Bmat = [A2, sparse(Ndof,Ndof);
    sparse(Ndof,Ndof), I_N];

nev  = mode_plot_number*2;
opts.tol = 1e-10;  opts.maxit = 3000;  opts.isreal = false;
%fprintf('Solving QEP for m  (sigma = %.4f)...\n', m_guess);
[Z, Mdiag] = eigs(Amat, Bmat, nev, m_guess, opts);
m_vals = diag(Mdiag);

%% ── Filter to physically meaningful branch and sort by proximity ──────────
% Keep modes with predominantly real m (COMSOL's M_float is real by construction
% for a lossless closed cavity); small imaginary parts here indicate solver noise.
ok = abs(imag(m_vals)) < 0.02*abs(real(m_vals));
m_vals = m_vals(ok);
Z      = Z(:, ok);
[~, si] = sort(abs(real(m_vals) - m_guess));
m_vals  = m_vals(si);
Z       = Z(:, si);

fprintf('\n%4s  %14s  %10s\n', 'Mode', 'M_float', 'nearest int');
for k = 1:min(6,numel(m_vals))
    fprintf('%4d  %14.6f  %10d\n', k, real(m_vals(k)), round(real(m_vals(k))));
end

%% ── Extract the mode shape for the branch closest to m_guess ──────────────
% z = [m*x; x]  ->  the field DOF vector is the SECOND half
x_field = Z(Ndof+1:end, 1);
M_float = real(m_vals(1));
fprintf('\nAt lambda0 = %.2f um:  M_float = %.4f  (nearest integer m = %d)\n', ...
    lambda0*1e6, M_float, round(M_float));

% m_guess    = n_eff_est * k0 * Rc
n_effs = m_vals./(k0*R_central);

nplot = min(mode_plot_number, size(Z,2));
Gamma = [];
lp = [1 2; 2 3; 3 1];   % local edge node pairs
for k = 1:nplot
    % z = [m*x ; x]  ->  field DOF vector is the SECOND half
    field_k = Z(Ndof+1:end, k);
    m_k     = real(m_vals(k));

    % Effective index and frequency for this branch, for the plot title
    n_eff_k = m_k / (k0 *R_central);
    freq_k  = c0 /lambda0;   % fixed target frequency (all branches share it)

    % % Normalize for consistent plot scaling (arbitrary eigenvector norm)
    % field_k = field_k / max(abs(field_k));

    if compute_fields==1
        [Er, Ep, Ez, Gamma] = reconstruct(field_k, rho_all, z_all, elements, ...
            eps_elem, Ndof_e, elemEdges, edgeSigns, lp, n_ring, n_clad);

        [Hr, Hp, Hz] = compute_Hfield(field_k, rho_all, z_all, elements, ...
            eps_elem, Ndof_e, elemEdges, edgeSigns, lp, n_ring, n_clad, m_k, lambda0);

        %% ── Power normalization ────────────────────────────────────────────────
        % Compute unnormalized guided power P0 = (1/2) Re{ integral (Er*Hp* - Ep*Hr*) dA }
        P0 = 0;
        for e = 1:numElem
            nd  = elements(e,:);
            re  = rho_all(nd)*1e6;  % convert um to meter
            ze = z_all(nd)*1e6;      % convert um to meter
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

        if plot_E_fields ==1
            plot_Emode(Er, Ep, Ez, Gamma, rho_all*1e6, z_all*1e6, elements, ...
                rho_inner*1e6, rho_outer*1e6, ring_height*1e6, k, freq_k, n_eff_k, ...
                save_figures, filename, k);
        end
        if plot_H_fields ==1
            plot_Hmode(Hr, Hp, Hz, Gamma, rho_all*1e6, z_all*1e6, elements, ...
                rho_inner*1e6, rho_outer*1e6, ring_height*1e6, k, freq_k, n_eff_k,...
                save_figures, [filename 'H'],k);
        end
    end

    fprintf('Mode %d:  M_float = %.4f,  n_eff = %.4f,  Gamma = %.4f\n', ...
        k, m_k, n_eff_k, Gamma);
end