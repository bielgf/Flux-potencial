function [Cl_y_w, Cl_y_h, Cdi_y_w, Cdp_y_w, Cd_y_w, Cdi_y_h, Cdp_y_h, Cd_y_h, CM_loc] = computeWingCanardAerodynamics(gamma_w,gamma_h,alpha,cwi05,chi05,Qinf_mod,Cl0w,Cl0h,Claw,Clah,thetai05,i_h,rho,dy_w,dy_h,l_h,Cm14_w,Cm14_h)

% --- Wing ---
Cl_y_w = (2 * gamma_w) ./ (Qinf_mod * cwi05);
alpha_i_w = (alpha + thetai05) - (Cl_y_w - Cl0w)/Claw; 
Cdi_y_w = Cl_y_w .* alpha_i_w;
Cdp_y_w = 0.0183 * Cl_y_w.^2 - 0.0302 * Cl_y_w + 0.0187; 
Cd_y_w = Cdi_y_w + Cdp_y_w;

% --- Canard ---
Cl_y_h = (2 * gamma_h) ./ (Qinf_mod * chi05);
alpha_i_h = (alpha + i_h) - (Cl_y_h - Cl0h)/Clah;
Cdi_y_h = Cl_y_h .* alpha_i_h;
Cdp_y_h = 0.0052 * Cl_y_h.^2 + 0.0071;
Cd_y_h = Cdi_y_h + Cdp_y_h;

% Integration of global forces and moments 
L_w = sum(rho * Qinf_mod * gamma_w .* dy_w);
L_h = sum(rho * Qinf_mod * gamma_h .* dy_h);

M14_w = sum(0.5 * rho * Qinf_mod^2 * cwi05.^2 * Cm14_w .* dy_w);
M14_h = sum(0.5 * rho * Qinf_mod^2 * chi05.^2 * Cm14_h .* dy_h);

% CM location (M_CM = 0)
if Claw ~= 0 && Clah ~= 0
    CM_loc = (M14_w + M14_h + L_h * l_h) / (L_w + L_h);
    fprintf('--- RESULTS (PART 2, SECTION 2) ---\n');
    fprintf('CM required position is: %.4f meters\n\n', CM_loc);
else
    CM_loc = 0;
end

end 