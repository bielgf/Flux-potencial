function polarAerodynamicPlot(Nw,Nh,P_w,P_h,P_w_mid,P_h_mid,cwi05,chi05,Qinf_mod,Cl0w,Cl0h,Claw,Clah,thetai05,i_h,ur,k_inf,rho,dy_w,dy_h,S_w,S_v)

alpha_range = 0:1:6;
CL_total = zeros(size(alpha_range));
CD_total = zeros(size(alpha_range));

for i = 1:length(alpha_range)
    alpha_i = deg2rad(alpha_range(i));
    
    [~, g_w, g_h] = computeWingCanardConfig(Nw,Nh,P_w,P_h,P_w_mid,P_h_mid,...
                    cwi05,chi05,Qinf_mod,Cl0w,Cl0h,Claw,Clah,thetai05,i_h,alpha_i,ur,k_inf);
                
    Cl_w = (2 * g_w) ./ (Qinf_mod * cwi05);
    al_i_w = (alpha_i + thetai05) - (Cl_w - Cl0w)/Claw;
    
    Cl_h = (2 * g_h) ./ (Qinf_mod * chi05);
    al_i_h = (alpha_i + i_h) - (Cl_h - Cl0h)/Clah;
    
    L_w_total = sum(rho * Qinf_mod * g_w .* dy_w);
    L_h_total = sum(rho * Qinf_mod * g_h .* dy_h);
    
    D_w_total = sum(rho * Qinf_mod * g_w .* al_i_w .* dy_w) + sum(0.5*rho*Qinf_mod^2*cwi05.*(0.0183*Cl_w.^2 - 0.0302*Cl_w + 0.0187).*dy_w);
    D_h_total = sum(rho * Qinf_mod * g_h .* al_i_h .* dy_h) + sum(0.5*rho*Qinf_mod^2*chi05.*(0.0052*Cl_h.^2 + 0.0071).*dy_h);
    D_vtp = 0.5 * rho * Qinf_mod^2 * S_v * 0.0062; % Resistencia simétrica VTP
    
    CL_total(i) = (L_w_total + L_h_total) / (0.5 * rho * Qinf_mod^2 * S_w);
    CD_total(i) = (D_w_total + D_h_total + D_vtp) / (0.5 * rho * Qinf_mod^2 * S_w);
end

% POLAR PLOT
figure;
plot(CD_total, CL_total, 'o-', 'LineWidth', 1.5);
xlabel('C_D (Global Drag Coefficient)');
ylabel('C_L (Global Lift Coefficient)');
title('Aerodynamic Polar Curve of the Glider (\alpha from 0° to 6°])');
grid on;

end 