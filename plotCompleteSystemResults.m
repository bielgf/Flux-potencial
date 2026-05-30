function plotCompleteSystemResults( ...
    P_w_mid, P_h_mid, b, ...
    gamma_w,     gamma_h,     gamma_w_iso,     gamma_h_iso, ...
    Cl_y_w,      Cl_y_h,      Cl_y_w_iso,      Cl_y_h_iso, ...
    Cdi_y_w,     Cdi_y_h,     Cdi_y_w_iso,     Cdi_y_h_iso, ...
    Cdp_y_w,     Cdp_y_h,     Cdp_y_w_iso,     Cdp_y_h_iso, ...
    Cd_y_w,      Cd_y_h,      Cd_y_w_iso,      Cd_y_h_iso)

    % Normalised spanwise coordinate
    eta_w = 2*P_w_mid(:,2)/b;
    eta_h = 2*P_h_mid(:,2)/b;
     
    % ------------------------------------------------------------------
    % Figure 1: Gamma and Cl
    % ------------------------------------------------------------------
    figure('Units', 'normalized', 'Position', [0.1, 0.2, 0.8, 0.45]);
    t = tiledlayout(1,2,'TileSpacing','compact');
     
    % --- TILE 1: Gamma ---
    ax1 = nexttile;
    hold on; box on;
    plot(eta_w, gamma_w,     'b');
    plot(eta_h, gamma_h,     'r');
    plot(eta_w, gamma_w_iso, 'b--');
    plot(eta_h, gamma_h_iso, 'r--');
    xlabel('2y/b');
    ylabel('Circulation, \Gamma [m^2/s]');
    title('Circulation distribution (\Gamma)');
    grid on;
    set(gca, 'GridAlpha', 0.15);
     
    % --- TILE 2: Cl ---
    ax2 = nexttile;
    hold on; box on;
    plot(eta_w, Cl_y_w,     'b');
    plot(eta_h, Cl_y_h,     'r');
    plot(eta_w, Cl_y_w_iso, 'b--');
    plot(eta_h, Cl_y_h_iso, 'r--');
    xlabel('2y/b');
    ylabel('Local lift coefficient, C_l');
    title('Lift distribution (C_l)');
    grid on;
    set(gca, 'GridAlpha', 0.15);
     
    lgd = legend(ax2, ...
        'Wing (Complete System)', ...
        'Canard (Complete System)', ...
        'Wing (Isolated)', ...
        'Canard (Isolated)');
    lgd.Layout.Tile = 'east';
     
    % ------------------------------------------------------------------
    % Figure 2: Cdi, Cdp and Cd
    % ------------------------------------------------------------------
    figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);
    t = tiledlayout(2,2,'TileSpacing','compact');
     
    % --- TILE 1: Cdi ---
    ax1 = nexttile;
    hold on; box on;
    plot(eta_w, Cdi_y_w,     'b');
    plot(eta_h, Cdi_y_h,     'r');
    plot(eta_w, Cdi_y_w_iso, 'b--');
    plot(eta_h, Cdi_y_h_iso, 'r--');
    xlabel('2y/b');
    ylabel('Local induced drag coefficient, C_{di}');
    title('Induced drag distribution (C_{di})');
    grid on;
    set(gca, 'GridAlpha', 0.15);
     
    % --- TILE 2: Cdp ---
    ax2 = nexttile;
    hold on; box on;
    plot(eta_w, Cdp_y_w,     'b');
    plot(eta_h, Cdp_y_h,     'r');
    plot(eta_w, Cdp_y_w_iso, 'b--');
    plot(eta_h, Cdp_y_h_iso, 'r--');
    xlabel('2y/b');
    ylabel('Local viscous drag coefficient, C_{dp}');
    title('Viscous drag distribution (C_{dp})');
    grid on;
    set(gca, 'GridAlpha', 0.15);
     
    % --- TILE 3: Cd ---
    ax3 = nexttile;
    hold on; box on;
    plot(eta_w, Cd_y_w,     'b');
    plot(eta_h, Cd_y_h,     'r');
    plot(eta_w, Cd_y_w_iso, 'b--');
    plot(eta_h, Cd_y_h_iso, 'r--');
    xlabel('2y/b');
    ylabel('Local drag coefficient, C_d');
    title('Drag distribution (C_d)');
    grid on;
    set(gca, 'GridAlpha', 0.15);
     
    lgd = legend(ax3, ...
        'Wing (Complete System)', ...
        'Canard (Complete System)', ...
        'Wing (Isolated)', ...
        'Canard (Isolated)');
    lgd.Layout.Tile = 4;
 
end