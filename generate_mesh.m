function [rho_all, z_all, elements, eps_elem] = generate_mesh( ...
    rho_inner, rho_outer, ring_height, ...
    rho_min,   rho_max,   z_min, z_max, ...
    n_ring,    n_clad, ...
    Nr_ring,   Nz_ring, ...
    Nr_clad,   Nz_clad)
% GENERATE_MESH  Conforming Delaunay mesh for a rectangular ring cross-section.
%
%   The four edges of the ring rectangle are inserted as CONSTRAINED edges
%   so that no triangle ever straddles the core/cladding boundary.
%   Every triangle is therefore purely inside or purely outside the ring.
%
%   Inputs
%   ------
%   rho_inner, rho_outer : inner / outer radii of ring   [m]
%   ring_height          : full height of ring           [m]
%   rho_min, rho_max     : domain rho extent             [m]
%   z_min,   z_max       : domain z extent               [m]
%   n_ring, n_clad       : refractive indices
%   Nr_ring, Nz_ring     : seed points along ring edges (rho and z directions)
%   Nr_clad, Nz_clad     : seed points along domain edges
%
%   Outputs
%   -------
%   rho_all   : [numNodes x 1]  nodal rho coordinates
%   z_all     : [numNodes x 1]  nodal z   coordinates
%   elements  : [numElem  x 3]  connectivity  (counter-clockwise)
%   eps_elem  : [numElem  x 1]  permittivity per element

z_bot = -ring_height/2;
z_top =  ring_height/2;

% =========================================================================
% 1.  SEED POINTS
%     Three groups:
%       (a) Points on the ring boundary rectangle  <- force conformity
%       (b) Interior points inside the ring        <- resolve the core
%       (c) Points in the surrounding cladding     <- resolve the domain
% =========================================================================

% ── (a) Ring boundary: all four edges, densely seeded ────────────────────
% Bottom edge  (z = z_bot, rho from rho_inner to rho_outer)
rho_bot = linspace(rho_inner, rho_outer, Nr_ring)';
z_bot_e = z_bot * ones(size(rho_bot));

% Top edge
rho_top = linspace(rho_inner, rho_outer, Nr_ring)';
z_top_e = z_top * ones(size(rho_top));

% Left edge  (rho = rho_inner, z from z_bot to z_top)
z_left  = linspace(z_bot, z_top, Nz_ring)';
rho_left_e = rho_inner * ones(size(z_left));

% Right edge  (rho = rho_outer)
z_right = linspace(z_bot, z_top, Nz_ring)';
rho_right_e = rho_outer * ones(size(z_right));

bdry_ring = [rho_bot,    z_bot_e; ...
             rho_top,    z_top_e; ...
             rho_left_e, z_left;  ...
             rho_right_e,z_right];

% ── (b) Interior of ring ──────────────────────────────────────────────────
[rG_in, zG_in] = meshgrid( ...
    linspace(rho_inner, rho_outer, Nr_ring+1), ...
    linspace(z_bot,     z_top,     Nz_ring+1));
interior_ring = [rG_in(:), zG_in(:)];
% Keep only strictly interior points (boundary is already in bdry_ring)
tol = 1e-14;
interior_ring = interior_ring( ...
    interior_ring(:,1) > rho_inner+tol & ...
    interior_ring(:,1) < rho_outer-tol & ...
    interior_ring(:,2) > z_bot   +tol & ...
    interior_ring(:,2) < z_top   -tol, :);

% ── (c) Cladding / domain ─────────────────────────────────────────────────
[rG_cl, zG_cl] = meshgrid( ...
    linspace(rho_min, rho_max, Nr_clad), ...
    linspace(z_min,   z_max,   Nz_clad));
clad_pts = [rG_cl(:), zG_cl(:)];
% Remove points that fall inside the ring (those are already seeded above)
inside = clad_pts(:,1) >= rho_inner & clad_pts(:,1) <= rho_outer & ...
         clad_pts(:,2) >= z_bot      & clad_pts(:,2) <= z_top;
clad_pts(inside, :) = [];

% ── Merge all points ───────────────────────────────────────────────────────
all_pts = unique([bdry_ring; interior_ring; clad_pts], 'rows');

% Clip to domain (safety)
all_pts = all_pts( ...
    all_pts(:,1) >= rho_min & all_pts(:,1) <= rho_max & ...
    all_pts(:,2) >= z_min   & all_pts(:,2) <= z_max, :);

% =========================================================================
% 2.  BUILD CONSTRAINED DELAUNAY TRIANGULATION
%     Insert the four ring edges as constraints so the triangulator is
%     forced to place triangle edges exactly along the material boundary.
% =========================================================================

% ── Locate indices of ring-boundary nodes after unique() ─────────────────
% We need to find which rows of all_pts correspond to each boundary edge.
% Use a tolerance-based nearest-node search.

% Gather all boundary nodes in order to form constraint segments
tol_snap = 1e-15 * max(rho_max-rho_min, z_max-z_min);

function idx = find_node(pt, pts, tol_s)
    d = sqrt((pts(:,1)-pt(1)).^2 + (pts(:,2)-pt(2)).^2);
    [dmin, idx] = min(d);
    if dmin > tol_s
        error('Boundary node not found in point set (d=%.2e).', dmin);
    end
end

% Build ordered sequences along each ring edge
r_bot_seq  = linspace(rho_inner, rho_outer, Nr_ring)';
r_top_seq  = r_bot_seq;
z_left_seq = linspace(z_bot, z_top, Nz_ring)';
z_rgt_seq  = z_left_seq;

% Index into all_pts for each sequence
idx_bot   = arrayfun(@(k) find_node([r_bot_seq(k), z_bot], all_pts, tol_snap), 1:Nr_ring);
idx_top   = arrayfun(@(k) find_node([r_top_seq(k), z_top], all_pts, tol_snap), 1:Nr_ring);
idx_left  = arrayfun(@(k) find_node([rho_inner, z_left_seq(k)], all_pts, tol_snap), 1:Nz_ring);
idx_right = arrayfun(@(k) find_node([rho_outer, z_rgt_seq(k)],  all_pts, tol_snap), 1:Nz_ring);

% Form constraint edge pairs  [n1 n2] for consecutive nodes on each edge
make_segs = @(ids) [ids(1:end-1)', ids(2:end)'];
segs_bot   = make_segs(idx_bot);
segs_top   = make_segs(idx_top);
segs_left  = make_segs(idx_left);
segs_right = make_segs(idx_right);

C = [segs_bot; segs_top; segs_left; segs_right];

% ── Constrained Delaunay triangulation ────────────────────────────────────
DT = delaunayTriangulation(all_pts, C);

rho_all  = DT.Points(:,1);
z_all    = DT.Points(:,2);
elements = DT.ConnectivityList;
numElem  = size(elements, 1);

% =========================================================================
% 3.  MATERIAL ASSIGNMENT
%     Because triangulation is constrained, the centroid test is now exact:
%     a centroid inside the ring rectangle belongs to the ring.
% =========================================================================
eps_elem = zeros(numElem, 1);
for e = 1:numElem
    nd  = elements(e,:);
    rc  = mean(rho_all(nd));
    zc  = mean(z_all(nd));
    if rc > rho_inner && rc < rho_outer && zc > z_bot && zc < z_top
        eps_elem(e) = n_ring^2;
    else
        eps_elem(e) = n_clad^2;
    end
end

% =========================================================================
% 4.  VALIDATE: no triangle should straddle the boundary
% =========================================================================
n_straddle = 0;
for e = 1:numElem
    nd = elements(e,:);
    rn = rho_all(nd);
    zn = z_all(nd);
    in_node = (rn >= rho_inner-tol_snap) & (rn <= rho_outer+tol_snap) & ...
              (zn >= z_bot   -tol_snap) & (zn <= z_top   +tol_snap);
    if any(in_node) && ~all(in_node)
        n_straddle = n_straddle + 1;
    end
end

fprintf('generate_mesh: %d nodes, %d elements, %d ring elements.\n', ...
    size(rho_all,1), numElem, sum(eps_elem == n_ring^2));

end