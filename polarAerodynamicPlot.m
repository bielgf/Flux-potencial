function polarAerodynamicPlot(Nw,Nh,P_w,P_h,P_w_mid,P_h_mid,cwi05,chi05,Qinf_mod,Cl0w,Cl0h,Claw,Clah,thetai05,i_h,ur,k_inf,rho,dy_w,dy_h,S_w,S_v)

alpha_range_deg = 0:0.5:6;
N_points = length(alpha_range_deg);

CL_total = zeros(1, N_points);
CD_total = zeros(1, N_points);

for i = 1:N_points
    alpha_i = deg2rad(alpha_range_deg(i));
    
    [~, g_w, g_h] = computeWingCanardConfig(Nw,Nh,P_w,P_h,P_w_mid,P_h_mid,...
                    cwi05,chi05,Qinf_mod,Cl0w,Cl0h,Claw,Clah,thetai05,i_h,alpha_i,ur,k_inf);
                
    Cl_w = (2 * g_w) ./ (Qinf_mod * cwi05);
    al_i_w = (alpha_i + thetai05) - (Cl_w - Cl0w)/Claw;
    
    Cl_h = (2 * g_h) ./ (Qinf_mod * chi05);
    al_i_h = (alpha_i + i_h) - (Cl_h - Cl0h)/Clah;
    
    L_w_total = sum(rho * Qinf_mod * g_w .* dy_w);
    L_h_total = sum(rho * Qinf_mod * g_h .* dy_h);
    
    D_w_ind  = sum(rho * Qinf_mod * g_w .* al_i_w .* dy_w);
    D_w_visc = sum(0.5*rho*Qinf_mod^2*cwi05.*(0.0183*Cl_w.^2 - 0.0302*Cl_w + 0.0187).*dy_w);
    
    D_h_ind  = sum(rho * Qinf_mod * g_h .* al_i_h .* dy_h);
    D_h_visc = sum(0.5*rho*Qinf_mod^2*chi05.*(0.0052*Cl_h.^2 + 0.0071).*dy_h);

    D_vtp    = 0.5 * rho * Qinf_mod^2 * S_v * 0.0062; 
    
    CL_total(i) = (L_w_total + L_h_total) / (0.5 * rho * Qinf_mod^2 * S_w);
    CD_total(i) = (D_w_ind + D_w_visc + D_h_ind + D_h_visc + D_vtp) / (0.5 * rho * Qinf_mod^2 * S_w);
end

% Aerodynamic Efficiency (E = L/D = CL/CD)
E_total = CL_total ./ CD_total;
[Emax_total, idx_max_total] = max(E_total);


% POLAR PLOT
figure('Units', 'normalized', 'Position', [0.15, 0.2, 0.7, 0.45], 'Color', 'w');

% --- (CL vs CD) ---
subplot(1,2,1); hold on; box on;
plot(CD_total, CL_total);
plot(CD_total(idx_max_total), CL_total(idx_max_total), 'ks', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

xlabel('Global Drag Coefficient, C_D');
ylabel('Global Lift Coefficient, C_L');
title('Aerodynamic Polar Curve of the Glider (C_L vs C_D)');
grid on; set(gca, 'GridAlpha', 0.15);
legend('Complete System (W+C+VTP)', 'E_{max}', 'Location', 'SouthEast');

% --- (L/D vs alpha) ---
subplot(1,2,2); hold on; box on;
plot(alpha_range_deg, E_total);
plot([0, alpha_range_deg(idx_max_total)], [Emax_total, Emax_total], 'k:', 'LineWidth', 1);
plot([alpha_range_deg(idx_max_total), alpha_range_deg(idx_max_total)], [0, Emax_total], 'k:', 'LineWidth', 1);
xline(alpha_range_deg(idx_max_total), '--', 'Color', 'k', 'LineWidth', 1.5, 'Label',...
    sprintf('E_{max} = %.2f at \\alpha = %.1fº\n', Emax_total, alpha_range_deg(idx_max_total)), 'LabelVerticalAlignment', 'middle');

xlabel('Angle of Attack, \alpha [º]');
ylabel('Efficiency, E = C_L / C_D'); ylim([0; 45]);
title('Aerodynamic Efficiency vs \alpha');
grid on; set(gca, 'GridAlpha', 0.15);

end
