
% =========================================================================
%  PLOT ONE MODE: |H_rho|, |H_phi|, |H_z| side by side
% =========================================================================
function plot_Hmode(Hr, Hp, Hz, Gamma, rho_um, z_um, elements, ...
    ri, ro, rh, modeNum, freq_val, neff_val,...
    save_figures, filename,k)

set(0,'defaultlinelinewidth',2)
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',18)

bx = [ri ri ro ro ri];
bz = [-rh rh rh -rh -rh]/2;

fields = {abs(Hr), abs(Hp), abs(Hz)};
labels = {'|{\it{H}}_\rho|', '|{\it{H}}_\phi|', '|{\it{H}}_z|'};
cmax   = max([max(abs(Hr)), max(abs(Hp)), max(abs(Hz))]);
if cmax == 0, cmax = 1; end

figure('Name', sprintf('Mode %d', modeNum), 'Position', [30 30 1380 430]);
for p = 1:3
    ax = subplot(1,3,p);
    trisurf(elements, rho_um, z_um, fields{p}, 'EdgeColor','none');
    view(2);  shading interp;  axis equal tight;  colorbar;
    %clim([0 cmax]);  
    colormap(ax, hot(256));
    yticklabels(strrep(yticklabels,'-','–'));
    xlabel('\rho [\mum]');  ylabel('z [\mum]');  title(labels{p});
    %hold on;  plot(bx, bz, 'c-', 'LineWidth', 1.8);
end
sgtitle(sprintf('Mode %d  |  f = %.4f THz  |  n_{eff} = %.4f  |  \\Gamma = %.3f', ...
    modeNum, freq_val/1e12, neff_val, Gamma));
%
if save_figures==1
    mfname = [filename '_Hmode' int2str(k)];
    print(gcf, mfname, '-dpng', '-r300');
end

xleft = min(rho_um)+(max(rho_um)-min(rho_um))/10;
yleft = max(z_um)*0.6;

figure('Name', sprintf('Mode %d', modeNum), 'Position', [530 30 500 900]);
for p = 1:3
    ax = subplot(3,1,p);
    trisurf(elements, rho_um, z_um, fields{p}, 'EdgeColor','none');
    view(2);  shading interp;  axis equal tight;  colorbar; 
    yticklabels(strrep(yticklabels,'-','–'));
    colormap(ax, hot(256));
    xlabel('\rho [\mum]');  ylabel('z [\mum]'); 
    text(xleft,yleft, labels{p},'color','w');
end

if save_figures==1
    mfname = [filename '_Hmode3x1_' int2str(k)];
    print(gcf, mfname, '-dpng', '-r300');
end

end

