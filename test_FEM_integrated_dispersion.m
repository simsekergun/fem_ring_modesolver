clear; close all; clc;

%% Parallel setup
% if your machine doesn't have 4 nodes or more, uncomment below
% replace the "parfor" in line 47 with "for"
parpool_no = 4;
p = gcp('nocreate');
if isempty(p)
    parpool(parpool_no);
elseif p.NumWorkers ~= parpool_no
    delete(gcp('nocreate'));
    parpool(parpool_no);
end

%% ── Constants ─────────────────────────────────────────────────────────────
c0  = 299792458;
mu0 = 4*pi*1e-7;
eps0 = 1/(mu0*c0^2);

%% ── Parameters ────────────────────────────────────────────────────────────
filename = 'test_Dint';
mesh_plot = 0;          % 1 means plot
mode_plot = 0;          % 1 means plot
save_figures = 0;       % 1 means save the figures

%% ── Geometry ──────────────────────────────────────────────────────────────
R_central   = 23e-6;
ring_width  = 890e-9;
ring_height = 670e-9;
material_ring = 'si3n4';
material_clad = 'sio2';

% Computational domain: padding on each side
pad_rho  = 2*ring_width;
pad_z    = 2*ring_height;

pump_freq = 283e12;
lambda_target = c0/pump_freq;              % pump wavelength
lambdas = (750:10:1500)*1e-9;

n_effs1 = zeros(length(lambdas),1);
n_effs2 = zeros(length(lambdas),1);

meshEPW = 30;  % elements per wavelength (higher == more accurate)

tic 
parfor counter = 1:length(lambdas)
    lambda0    = lambdas(counter);   
   
    [n_effs1(counter), n_effs2(counter)] = FEM_Ring_Solver_neff_interp(lambda0, material_ring, ...
        material_clad, R_central, ring_width, ring_height, meshEPW, pad_rho, pad_z, mode_plot);
end
toc

%%% Integrated Dispersion Calculations %%%
[muA, DintA, D1A, D2A, D3A, D4A] = anaylze_dispersion(R_central, lambdas, n_effs1.',pump_freq);
[muB, DintB, D1B, D2B, D3B, D4B] = anaylze_dispersion(R_central, lambdas, n_effs2.',pump_freq);

save('results_FEM_test890_670.mat','muA','muB','lambdas','DintA','DintB','n_effs1','n_effs2');

% plot integrated dispersion results
fig = figure(2); clf;
fig.Position=[100 100 800 600];
subplot(211); 
plot(lambdas*1e6, n_effs1,lambdas*1e6, n_effs2)
xlim([min(lambdas*1e6) max(lambdas*1e6)]);
legend('Mode 1','Mode 2');
xlabel('Wavelength (\mum)');
ylabel('Effective Index');
grid on;

subplot(212); 
plot(muA, DintA/1e9,muB, DintB/1e9);
legend('Mode 1','Mode 2','Location','SouthWest');
xlim([min([muA muB]) max([muA muB])])
xlabel('Mode Index');
ylabel('{\it{D}}_{int} (GHz)');
grid on;
print -dpng figure_Dint_FEM_890_670
