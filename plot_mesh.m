function plot_mesh(rho_all, z_all, elements, eps_elem, rho_inner, rho_outer, ring_height,save_figures, filename)
% PLOT_MESH  Visualize the FEM mesh with triangles colored by permittivity.
%
%   plot_mesh(rho_all, z_all, elements, eps_elem, rho_inner, rho_outer, ring_height)
%
%   Inputs
%   ------
%   rho_all     : [numNodes x 1]  nodal rho coordinates  [m]
%   z_all       : [numNodes x 1]  nodal z   coordinates  [m]
%   elements    : [numElem  x 3]  node connectivity
%   eps_elem    : [numElem  x 1]  relative permittivity per element
%   rho_inner   : inner radius of ring  [m]
%   rho_outer   : outer radius of ring  [m]
%   ring_height : height of ring cross-section  [m]
set(0,'defaultlinelinewidth',3)
set(0,'DefaultAxesFontSize',24)
set(0,'DefaultTextFontSize',24)

rho_um = rho_all * 1e6;   % convert to microns for display
z_um   = z_all   * 1e6;

numElem    = size(elements, 1);
eps_unique = unique(eps_elem);
nMat       = numel(eps_unique);

% ── Assign a color to each unique permittivity value ─────────────────────
% Use a perceptually distinct colormap: blue for cladding, orange for ring
cmap_mat = [0.25 0.55 0.85;   % blue   -> lower eps (cladding)
            0.95 0.50 0.10;   % orange -> higher eps (ring)
            0.35 0.75 0.35;   % green  -> third material if present
            0.80 0.25 0.25];  % red    -> fourth material if present
cmap_mat = cmap_mat(1:min(nMat, size(cmap_mat,1)), :);

% Map each element to its material color
faceColors = zeros(numElem, 3);
for k = 1:nMat
    idx = eps_elem == eps_unique(k);
    faceColors(idx, :) = repmat(cmap_mat(k,:), sum(idx), 1);
end

% ── Draw triangles ────────────────────────────────────────────────────────
figure('Name','FEM Mesh — Permittivity Distribution', ...
       'Position', [100 100 820 680]);
ax = axes;
hold(ax, 'on');

% patch() accepts per-face colors directly
xv = reshape(rho_um(elements'), 3, numElem);   % 3 x numElem
yv = reshape(z_um  (elements'), 3, numElem);

patch(ax, xv, yv, zeros(3,numElem), ...        % z-data unused in 2-D view
      'FaceVertexCData', faceColors, ...
      'FaceColor',       'flat', ...
      'EdgeColor',       [0.95 0.95 0.95], ...
      'LineWidth',        0.3);

% ── Ring outline ──────────────────────────────────────────────────────────
ri = rho_inner * 1e6;
ro = rho_outer * 1e6;
rh = ring_height / 2 * 1e6;
bx = [ri ri ro ro ri];
bz = [-rh  rh  rh -rh -rh];
plot(ax, bx, bz, 'w-', 'LineWidth', 1.8);

% ── Colorbar / legend ─────────────────────────────────────────────────────
% Build a manual legend using invisible scatter points
hLeg = gobjects(nMat, 1);
for k = 1:nMat
    hLeg(k) = scatter(ax, nan, nan, 60, cmap_mat(k,:), 'filled');
end
legend(ax, hLeg, arrayfun(@(e) sprintf('\\epsilon_r = %.4g', e), ...
       eps_unique, 'UniformOutput', false), ...
       'Location', 'NorthEast');

% ── Labels ────────────────────────────────────────────────────────────────
axis(ax, 'equal', 'tight');
xlabel(ax, '\rho  [\mum]');
ylabel(ax, 'z  [\mum]');
%title(ax,  'FEM Mesh — Permittivity Distribution', 'FontSize', 12);

% Annotation: node and element count
numNodes = size(rho_all, 1);
numEdges_est = round(1.5 * numElem);   % Euler estimate
annotation('textbox', [0.13 0.01 0.5 0.04], ...
    'String', sprintf('Nodes: %d   |   Elements: %d   |   Materials: %d', ...
                      numNodes, numElem, nMat), ...
    'EdgeColor', 'none', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);

view(ax, 2);
box(ax, 'on');
hold(ax, 'off');
yticklabels(strrep(yticklabels,'-','–'));

if save_figures==1
    mfname = [filename '_mesh'];
    print(gcf, mfname, '-dpng', '-r300');
end

end