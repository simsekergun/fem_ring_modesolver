
% =========================================================================
%  PLOT ONE MODE: |E_rho|, |E_phi|, |E_z| side by side
% =========================================================================
function plot_mode(Er, Ep, Ez, Gamma, rho_all, z_all, elements, ...
    rho_inner, rho_outer, ring_height, modeNum, freq_val, neff_val,...
    save_figures, filename,k)

set(0,'defaultlinelinewidth',2)
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',18)

rho_um = rho_all*1e6;  z_um = z_all*1e6;
ri = rho_inner*1e6;    ro = rho_outer*1e6;
rh = ring_height/2*1e6;
bx = [ri ri ro ro ri];
bz = [-rh rh rh -rh -rh];

fields = {abs(Er), abs(Ep), abs(Ez)};
labels = {'|{\it{E}}_\rho|', '|{\it{E}}_\phi|', '|{\it{E}}_z|'};
cmax   = max([max(abs(Er)), max(abs(Ep)), max(abs(Ez))]);
if cmax == 0, cmax = 1; end

figure('Name', sprintf('Mode %d', modeNum), 'Position', [30 30 1380 430]);
for p = 1:3
    ax = subplot(1,3,p);
    trisurf(elements, rho_um, z_um, fields{p}, 'EdgeColor','none');
    view(2);  shading interp;  axis equal tight;  colorbar;
    colormap(ax, hot(256));
    xlabel('\rho [\mum]');  ylabel('z [\mum]');  title(labels{p});
    hold on;  plot(bx, bz, 'c-', 'LineWidth', 1.8);
end
sgtitle(sprintf('Mode %d  |  f = %.4f THz  |  n_{eff} = %.4f  |  \\Gamma = %.3f', ...
    modeNum, freq_val/1e12, neff_val, Gamma));
%
if save_figures==1
    mfname = [filename '_mode' int2str(k)];
    print(gcf, mfname, '-dpng', '-r300');
end



%% un-comment below, if you want to plot individually
% for p = 1:3    
%     figure('Name', sprintf('Mode %d %d', modeNum, p), 'Position', [300 300 600 400]);
%     trisurf(elements, rho_um, z_um, fields{p}, 'EdgeColor','none');
%     view(2);  shading interp;  axis equal tight;  colorbar;
%     %clim([0 cmax]);  
%     colormap(hot(256));
%     xlabel('\rho [\mum]');  ylabel('z [\mum]');      
%     text(min(rho_um)*1.02, max(z_um)*0.8, labels{p},'Color','w');
%     hold on;  plot(bx, bz, 'c-', 'LineWidth', 1.8);
%     yticklabels(strrep(yticklabels,'-','–'));
%     if save_figures==1
%         mfname = ['./figures/' filename '_mode' int2str(k) int2str(p)];
%         print(gcf, mfname, '-dpng', '-r300');
%     end
% end

end

