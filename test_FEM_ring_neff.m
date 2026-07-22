% =========================================================================
% AXISYMMETRIC DIELECTRIC RING RESONATOR — FEM EIGENMODE SOLVER
% =========================================================================
clear; close all; clc;

%% ── input parameters
filename = 'test_RectCS';
mesh_plot = 1;          % 1 means plot
mode_plot = 1;          % 1 means plot
number_of_modes = 2;
save_figures = 1;       %  1 means save the figures

lambda0    = 1060e-9;       % wavelength

%% materials
ring_material = 'si3n4';
clad_material = 'sio2';

%% Geometry 
R_central   = 23e-6;
ring_width  = 890e-9;
ring_height = 670e-9;
%% Computational domain: padding on each side
pad_rho  = 2*ring_width;
pad_z    = 2*ring_height;

%% mesh quality (higher EPW = higher accuracy)
% EPW of 40 requires >1 GB of RAM
% use a smaller value if your machine has a limited memory
meshEPW = 20;   % elements per wavelength

%% call the main solver (assuming a rectangular cross-section)
ring_fem_solver;

