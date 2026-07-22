clear;
% all dimensions in micrometer
%% ── Constants ─────────────────────────────────────────────────────────────
c0  = 299792458;
mu0 = 4*pi*1e-7;
eps0 = 1/(mu0*c0^2);

%% materials
core_material = 'si3n4';
clad_material = 'sio2';

%% Target: fix wavelength, solve for m 
lambda0    = 1.06e-6;

%% ── Geometry ──────────────────────────────────────────────────────────────
R_central   = 23e-6;
ring_width  = 0.890e-6;
ring_height = 0.670e-6;
pad_rho  = 2*ring_width;
pad_z    = 2*ring_height;

meshEPW = 20;

filename     = 'test_solve_for_m';      % to save the outputs
mode_plot_number = 4;      % how many top modes to plot
compute_fields = 1;         % if you want to output the fields
% visualization settings
plot_E_fields = 1;                 % if you also want to plot the E fields
plot_H_fields = 1;                 % if you also want to plot the H fields
mesh_plot  = 1;                   % set to 1 if you want to visualie the mesh
save_figures = 1;                 % set to 1 to save PNGs
save_results = 1;

%% run the solver
ring_rcs_m_solver