%------------------------------------------------------------------%
%%%%% Flux Potencial. Aerodinàmica Numèrica en Perfils i Ales %%%%%%
%%%%%%%%%%%%%%% Judith Bailen, Eulàlia Caballol %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Biel González, Júlia Soliva  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% AMVO - Course 2025-2026 %%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------%


% ---------------------------------------------------------------- %
% ---------- PART 1: THE CONSTANT STRENGTH VORTEX METHOD --------- %
% ---------- APPLIED TO SIMPLE AND TWO-ELEMENT AIRFOILS ---------- %
% ---------------------------------------------------------------- %

clear;clc;close all

%% ------------------ HQ300 AIRFOIL STUDY ------------------- %%

alpha = deg2rad(0:2:8);
Ndiv  = [16, 32, 64, 128, 256, 512];

nA = length(alpha);
nN = length(Ndiv);

CL_table   = zeros(nA,nN);
CM14_table = zeros(nA,nN);
L_table    = zeros(nA,nN);
M14_table  = zeros(nA,nN);

CM0      = zeros(1,nA);  
AC_table = zeros(1,nN);

rho     = 1.225;
gam     = 1.4;
Qinfmod = 1;
Minf    = 0;


for jj = 1:nN
    for ii = 1:nA
        Qinf  = Qinfmod*[cos(alpha(ii)),sin(alpha(ii))];
        geom  = computeGeometry(Ndiv(jj));
        gamma = computeCSV(Ndiv(jj),geom,Qinf);
        [CL,L,CM1_4,M1_4,CM0(ii),cp,~] = computeAerodynamics(Ndiv(jj),geom,gamma,Qinfmod,rho,alpha(ii),Minf);
        CL_table(ii,jj)   = CL;
        CM14_table(ii,jj) = CM1_4;
        L_table(ii,jj)    = L;
        M14_table(ii,jj)  = M1_4;
        
        if Ndiv(jj) == 32  && alpha(ii) == deg2rad(8); plotGammaDistribution(Ndiv(jj),geom.X,geom.Nc,gamma); end
        if Ndiv(jj) == 256 && alpha(ii) == deg2rad(8); plotCpChordDistribution(geom.X,cp,alpha(ii),Ndiv(jj));
                                                       plotCpDistribution(geom.X,geom.Nc,geom.Xc,cp,alpha(ii),Ndiv(jj),CL); end
    end
    dCM0_da = polyfit(alpha(:),CM0(:),1);
    dCL_da  = polyfit(alpha(:),CL_table(:,jj),1);
    AC_table(jj) = -dCM0_da(1)/dCL_da(1);
end

plotCLandCM14vsalpha(alpha,CL_table,CM14_table);
plotConvergence(Ndiv,CL_table,alpha);


Mcr = zeros(1,nA-2);

for ii = 1:nA-2
    Qinf  = Qinfmod*[cos(alpha(ii)),sin(alpha(ii))];
    geom  = computeGeometry(512);
    gamma = computeCSV(512,geom,Qinf);
    [~,~,~,~,~,cp,~] = computeAerodynamics(512,geom,gamma,Qinfmod,rho,alpha(ii),Minf);
    Mcr(ii) = computeCriticalMachNumber(cp,gam);
end


Mcr_val      = [0, Mcr(end) - 0.15, Mcr(end) - 0.10, Mcr(end) - 0.05, Mcr(end)]; % The first being 0 to compare with the results obtained previously
CL_Mcr_table = zeros(1,length(Mcr_val));

for ii = 1:length(Mcr_val)
    Qinf  = Qinfmod*[cos(alpha(3)),sin(alpha(3))];
    geom  = computeGeometry(512);
    gamma = computeCSV(512,geom,Qinf);
    [~,~,~,~,~,~,CL_comp] = computeAerodynamics(512,geom,gamma,Qinfmod,rho,alpha(3),Mcr_val(ii));
    CL_Mcr_table(ii) = CL_comp;
end

%% ------------------ TWO NACA 0012 AIRFOILS TANDEM ------------------ %%
clear;clc;

% Disclaimer: It is important for this code to mantain the same number of
% Ndiv_main and Ndiv_second, as well as the same number of alpha and
% delta_e combinations, for how the loops are defined.

Ndiv_main   = [16, 32, 64, 128, 256, 512];
Ndiv_second = [16, 32, 64, 128, 256, 512];
Ntotal      = Ndiv_main + Ndiv_second;
subsection  = 5;  % 4 or 5

switch subsection
    case 4
        alpha   = deg2rad(0:2:8);
        delta_e = [0;0;0;0;0];
    case 5
        alpha   = deg2rad([4;4;4;4;4]);
        delta_e = 0:4:16;
    otherwise
        disp('Incorrect exercise subsection')
end

nA = length(alpha);
nN = length(Ndiv_main);

CL_table   = zeros(nA,nN);
CM14_table = zeros(nA,nN);
L_table    = zeros(nA,nN);
M14_table  = zeros(nA,nN);

rho     = 1.225;
gam     = 1.4;
Qinfmod = 1;
Minf    = 0;

for jj = 1:nN
    for ii = 1:nA
        Qinf  = Qinfmod*[cos(alpha(ii)),sin(alpha(ii))];
        geom  = computeGeometry2airfoil(Ndiv_main(jj), Ndiv_second(jj), delta_e(ii));

        gamma = computeCSV2(Ndiv_main(jj),Ndiv_second(jj),geom,Qinf);
        [CL,L,CM1_4,M1_4,CM0,cp,~] = computeAerodynamics2(Ndiv_main(jj),Ndiv_second(jj),geom,gamma,Qinfmod,rho,alpha(ii),Minf);
        CL_table(ii,jj)   = CL;
        CM14_table(ii,jj) = CM1_4;
        L_table(ii,jj)    = L;
        M14_table(ii,jj)  = M1_4;     
    end

    if Ndiv_main(jj) == 256 && subsection == 4; plotCpChordDistribution2(geom.X,cp,alpha(ii),delta_e(ii),Ndiv_main(jj),Ndiv_second(jj));
                                                plotCpDistribution2(geom.X,geom.Nc,geom.Xc,cp,alpha(ii),delta_e(ii),Ndiv_main(jj),Ndiv_second(jj)); end
    if Ndiv_main(jj) == 256 && subsection == 5; plotCpChordDistribution2(geom.X,cp,alpha(ii),delta_e(ii),Ndiv_main(jj),Ndiv_second(jj));
                                                plotCpDistribution2(geom.X,geom.Nc,geom.Xc,cp,alpha(ii),delta_e(ii),Ndiv_main(jj),Ndiv_second(jj)); end 
end

% PLOTS
if subsection == 4; plotCLandCM14vsalpha(alpha,CL_table,CM14_table); end
if subsection == 5; plotCLandCM14vsdelta(delta_e,CL_table,CM14_table); end

%% ------------------------------------------------------------------- %%
% --------------- PART 2: PRANDTL-S LIFTING LINE MODEL ---------------- %
% ---------- APPLIED TO COMPOUND WINGS OF LARGE ASPECT RATIO ---------- %
% --------------------------------------------------------------------- %

clc; clear;

% ---------------------- INPUT DATA ---------------------- %

% Geometric Data
b      = 15;
b_h    = 3;
c_r    = 0.95;
c_rh   = 0.5;
c_t    = 0.55;
c_th   = 0.3;
l_h    = -4;
l_v    = 1.2;
Sv     = 1.5;
i_w    = 0;
i_h    = deg2rad(3);
Sw     = (c_r + c_t)/2*b;
lambda = c_t/c_r;
c_bar  = (2/3)*c_r*(1 + lambda + lambda^2)/(1 + lambda); % MAC mean aerodynamic chord [m]

% Numerical Data
Nw = 512;
Nh = Nw/4;

% Aerodynamic Data
rho       = 1.225;
alpha     = deg2rad(4);
Claw      = 6.935845305522705;  
Cl0w      = 0.623151006958811;
i_inf     = [cos(alpha) , 0 , sin(alpha)];
ur        = -i_inf;
k_inf     = [-sin(alpha) , 0 , cos(alpha)];
Qinf_mod  = norm(i_inf);

Cd_w = @(Clw) 0.0183*Clw.^2 - 0.0302*Clw + 0.0187;
Cd_h = @(Clh) 0.0052*Clh^2 + 0.0071;
Cd_v = 0.0062;

%---------- GEOMETRIC DISCRETIZATION OF THE W-C ----------%

dy_w = b/Nw;
dy_h = b_h/Nh;
yP_w = -b/2:dy_w:b/2;
yP_h = -b_h/2:dy_h:b_h/2;

P_w     = [zeros(size(yP_w)) ; yP_w ; zeros(size(yP_w))]';
P_h     = [l_h*ones(size(yP_h)) ; yP_h ; zeros(size(yP_h))]';
P_w_mid = [zeros(1,size(yP_w,2)-1) ; (yP_w(1:end-1)+yP_w(2:end))/2 ; zeros(1,size(yP_w,2)-1)]';
P_h_mid = [l_h*ones(1,size(yP_h,2)-1) ; (yP_h(1:end-1)+yP_h(2:end))/2 ; zeros(1,size(yP_h,2)-1)]';

cwi05     = c_r + (c_t - c_r)*(2*abs(P_w_mid(:,2))/b);
chi05     = c_rh + (c_th - c_rh)*(2*abs(P_h_mid(:,2))/b_h);

%% ------------------ STUDY OF THE WING ISOLATED - HQ300 ------------------ %

twist_val_deg = [0:-0.25:-3 -3.01:-0.01:-3.99 -4:-0.25:-8];
twist_val     = deg2rad(twist_val_deg);
n_twist       = length(twist_val);

CL     = zeros(1, n_twist);
CD     = zeros(1, n_twist);
CD_ind = zeros(1, n_twist);
L      = zeros(1, n_twist);
D      = zeros(1, n_twist);

Cl_vec        = zeros(Nw,n_twist);
Cd_vec        = zeros(Nw,n_twist);
alpha_ind_vec = zeros(Nw, n_twist);
Cd_visc_vec   = zeros(Nw, n_twist);
Cd_ind_vec_1  = zeros(Nw, n_twist);
Cd_visc_vec_1 = zeros(Nw, n_twist);

for ii = 1:n_twist

    theta_mid = twist_val(ii)*(2*abs(P_w_mid(:,2))/b);
    gamma     = computeWingConfig(Nw,P_w,P_w_mid,ur,cwi05,k_inf,Claw,Qinf_mod,Cl0w,alpha,theta_mid);

    Cl_vec(:,ii) = (2*gamma)./(cwi05.*Qinf_mod);
    Cd_visc_vec  = Cd_w(Cl_vec(:,ii));
    alpha_ind    = (Cl_vec(:,ii) - Cl0w)/Claw - alpha - theta_mid;
    Cd_ind_vec   = -2*gamma.*alpha_ind./(Qinf_mod.*cwi05);
    
    Cd_vec(:,ii)        = Cd_visc_vec + Cd_ind_vec;
    Cd_ind_vec_1(:,ii)  = Cd_ind_vec;
    alpha_ind_vec(:,ii) = alpha_ind;          
    Cd_visc_vec_1(:,ii) = Cd_visc_vec; 
    
    CL(1,ii)     = sum(Cl_vec(:,ii).*cwi05.*(P_w(2:end,2) - P_w(1:end-1,2)))/Sw;
    CD(1,ii)     = sum(Cd_vec(:,ii).*cwi05.*(P_w(2:end,2) - P_w(1:end-1,2)))/Sw;
    CD_ind(1,ii) = sum(Cd_ind_vec.*cwi05.*(P_w(2:end,2) - P_w(1:end-1,2)))/Sw;
    L(1,ii)      = Qinf_mod*Sw*CL(1,ii);
    D(1,ii)      = Qinf_mod*Sw*CD(1,ii);

end

% Post - process
[~, idx_opt] = min(CD_ind);
theta_cd_min_deg = rad2deg(twist_val(idx_opt));
fprintf('--- RESULTS (PART 2, SECTION 1) ---\n');
fprintf('Twist que minimitza CDind: theta_t = %+.2f°\n', theta_cd_min_deg);

y_mid = (P_w(1:end-1,2) + P_w(2:end,2))/2;
eta   = y_mid/(b/2);   % coordenada normalitzada

theta_max_L_D = 0;
valor_max     = 0;

for jj = 1:n_twist
    valor = CL(jj)/CD(jj);
    if valor > valor_max
        valor_max = valor;
        theta_max_L_D = rad2deg(twist_val(jj));
    end
end

fprintf('Twist òptim (maximització L/D): theta_t = %+.2f°\n\n', theta_max_L_D);

plotIsolatedWingResults(twist_val_deg,twist_val,eta,Cl_vec,Cd_ind_vec_1,Cd_visc_vec_1,alpha_ind_vec,CL,CD_ind,CD,theta_cd_min_deg,theta_max_L_D,alpha);

%% ------------------ STUDY OF THE COMPLETE SYSTEM (W C VTP) ------------------ %%

Cl0h      = 0;                           % De la Part 1, Apartat 4
Clah      = (0.902768 - 0)/deg2rad(8);   % De la Part 1, Apartat 4
twist_tip = deg2rad(theta_max_L_D);
thetai05  = twist_tip*(2*abs(P_w_mid(:,2))/b);
Cm14_w    = -0.1376;                     % De la Part 1, Apartat 1: alpha 4
Cm14_h    = -0.0015;                     % De la Part 1, Apartat 5: alpha 4, delta 0


% 2. Spanwise distribution of aerodynamic coefficients and CM location
% (M_CM = 0) alpha = 4, delta = 0

alpha = deg2rad(4);
[gamma,gamma_w,gamma_h] = computeWingCanardConfig(Nw,Nh,P_w,P_h,P_w_mid,P_h_mid,cwi05,...
    chi05,Qinf_mod,Cl0w,Cl0h,Claw,Clah,thetai05,i_h,alpha,ur,k_inf);
[Cl_y_w, Cl_y_h, Cdi_y_w, Cdp_y_w, Cd_y_w, Cdi_y_h, Cdp_y_h, Cd_y_h, CM_loc] = computeWingCanardAerodynamics(gamma_w,gamma_h,...
    alpha,cwi05,chi05,Qinf_mod,Cl0w,Cl0h,Claw,Clah,thetai05,i_h,rho,dy_w,dy_h,l_h,Cm14_w,Cm14_h);

% Isolated wing
Cl0h_iso = 0; 
Clah_iso = 0;
[~, gamma_w_iso, ~] = computeWingCanardConfig(Nw,Nh,P_w,P_h,P_w_mid,P_h_mid,cwi05,chi05,...
    Qinf_mod,Cl0w,Cl0h_iso,Claw,Clah_iso,thetai05,i_h,alpha,ur,k_inf);
[Cl_y_w_iso, ~, Cdi_y_w_iso, Cdp_y_w_iso, Cd_y_w_iso, ~, ~, ~, ~] = computeWingCanardAerodynamics(gamma_w_iso,zeros(size(gamma_h)),...
    alpha,cwi05,chi05,Qinf_mod,Cl0w,Cl0h_iso,Claw,Clah_iso,thetai05,i_h,rho,dy_w,dy_h,l_h,Cm14_w,Cm14_h);

% Isolated canard
Cl0w_iso = 0; 
Claw_iso = 0;
[~, ~, gamma_h_iso] = computeWingCanardConfig(Nw,Nh,P_w,P_h,P_w_mid,P_h_mid,cwi05,chi05,...
    Qinf_mod,Cl0w_iso,Cl0h,Claw_iso,Clah,thetai05,i_h,alpha,ur,k_inf);
[~, Cl_y_h_iso, ~, ~, ~, Cdi_y_h_iso, Cdp_y_h_iso, Cd_y_h_iso, ~] = computeWingCanardAerodynamics(zeros(size(gamma_w)),gamma_h_iso,...
    alpha,cwi05,chi05,Qinf_mod,Cl0w_iso,Cl0h,Claw_iso,Clah,thetai05,i_h,rho,dy_w,dy_h,l_h,Cm14_w,Cm14_h);

plotCompleteSystemResults( ...
    P_w_mid, P_h_mid, b, ...
    gamma_w,     gamma_h,     gamma_w_iso,     gamma_h_iso, ...
    Cl_y_w,      Cl_y_h,      Cl_y_w_iso,      Cl_y_h_iso, ...
    Cdi_y_w,     Cdi_y_h,     Cdi_y_w_iso,     Cdi_y_h_iso, ...
    Cdp_y_w,     Cdp_y_h,     Cdp_y_w_iso,     Cdp_y_h_iso, ...
    Cd_y_w,      Cd_y_h,      Cd_y_w_iso,      Cd_y_h_iso);

% 3. Polar Aerodynamic Curve for delta = 0

polarAerodynamicPlot(Nw,Nh,P_w,P_h,P_w_mid,P_h_mid,cwi05,chi05,Qinf_mod,Cl0w,Cl0h,...
    Claw,Clah,thetai05,i_h,ur,k_inf,rho,dy_w,dy_h,Sw,Sv)


% 4. C_L and C_M for alpha = 4 and delta = 12

Cl_d12     = 1.3589;          % Part 1, Apartat 5: alpha 4, delta 12
Cm14_h_d12 = -0.1368;         % Part 1, Apartat 5: alpha 4, delta 12

alpha = deg2rad(4);
Cl0h_eff = Cl_d12 - Clah * alpha;  % elevator effect
[~, gamma_w_d, gamma_h_d] = computeWingCanardConfig(Nw,Nh,P_w,P_h,P_w_mid,P_h_mid,...
                             cwi05,chi05,Qinf_mod,Cl0w,Cl0h_eff,Claw,Clah,thetai05,i_h,alpha,ur,k_inf);

L_w_d         = sum(rho * Qinf_mod * gamma_w_d .* dy_w);
L_h_d         = sum(rho * Qinf_mod * gamma_h_d .* dy_h);
CL_global_d12 = (L_w_d + L_h_d) / (0.5 * rho * Qinf_mod^2 * Sw);

M14_w_d    = sum(0.5 * rho * Qinf_mod^2 * cwi05.^2 * Cm14_w .* dy_w);
M14_h_d    = sum(0.5 * rho * Qinf_mod^2 * chi05.^2 * Cm14_h_d12 .* dy_h); 
M_CM_total = M14_w_d - L_w_d * CM_loc + M14_h_d + L_h_d * (l_h - CM_loc);

Cm_global_d12 = M_CM_total / (0.5 * rho * Qinf_mod^2 * Sw * c_bar);

fprintf('--- RESULTS (PART 2, SECTION 4) ---\n');
fprintf('Global lift coefficient (C_L): %.4f\n', CL_global_d12);
fprintf('Pitching moment coefficient about CM (C_m,CM): %.4f\n', Cm_global_d12);