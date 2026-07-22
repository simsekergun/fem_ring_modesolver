% this subroutine calls get_ring_neff.m for different m values
% and calculates the effective indices of the first two resonant modes
% via a simple interpolation

function [n_eff1, n_eff2] = FEM_Ring_Solver_neff_interp(lambda0, material_ring, material_clad, ...
    R_central, ring_width, ring_height, meshEPW, pad_rho, pad_z, mode_plot)

k0         = 2*pi/lambda0;
n_ring     = get_refractive_index(material_ring,lambda0*1e6);
n_clad     = get_refractive_index(material_clad,lambda0*1e6);

numNodes = 2;
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
m_az = round(n_eff_est * k0 * R_central);

[rho_all, z_all, elements, eps_elem] = generate_mesh( ...
    rho_inner, rho_outer, ring_height, ...
    rho_min,   rho_max,   z_min, z_max, ...
    n_ring,    n_clad, ...
    Nr_ring,   Nz_ring, ...
    Nr_clad,   Nz_clad);


%% ── Edge table ────────────────────────────────────────────────────────────
[edges, elemEdges, edgeSigns] = build_edges(elements);
numEdges = size(edges,1);

Ndof_e = numEdges;
Ndof_p = numNodes;
Ndof   = Ndof_e + Ndof_p;

[n_effs1, n_effs2, lam_res1, lam_res2] = get_ring_neff(lambda0, ...
            R_central, ring_width, ring_height, n_ring, n_clad,rho_all, z_all, elements, eps_elem,...
            edges, elemEdges, edgeSigns,pad_rho,pad_z,mode_plot, m_az);        

m_az0 = m_az;
lam_res1b = lam_res1;

if lam_res1b <= lambda0*1e6
    while lam_res1b < lambda0*1e6

        disp('m value has been decreased by one...')
        m_az = m_az-1;
        [n_effs1b, n_effs2b, lam_res1b, lam_res2b] = get_ring_neff(lambda0, ...
            R_central, ring_width, ring_height, n_ring, n_clad,rho_all, z_all, elements, eps_elem,...
            edges, elemEdges, edgeSigns,pad_rho,pad_z,mode_plot, m_az);        
    end
else
    while lam_res1b >= lambda0*1e6
        disp('m value has been increased by one...')
        m_az = m_az+1;

        [n_effs1b, n_effs2b, lam_res1b, lam_res2b] = get_ring_neff(lambda0, ...
            R_central, ring_width, ring_height, n_ring, n_clad,rho_all, z_all, elements, eps_elem,...
            edges, elemEdges, edgeSigns,pad_rho,pad_z,mode_plot, m_az);        
    end
    
end


n_eff1 = n_effs1b + (lambda0*1e6 - lam_res1b) * (n_effs1 - n_effs1b) / (lam_res1 - lam_res1b);
n_eff2 = n_effs2b + (lambda0*1e6 - lam_res2b) * (n_effs2 - n_effs2b) / (lam_res2 - lam_res2b);