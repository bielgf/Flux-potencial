%------------------------------------------------------------------%
%%%%% Flux Potencial. Aerodinàmica Numèrica en Perfils i Ales %%%%%%
%%%%%%%%%%%%%%% Judith Bailen, Eulàlia Caballol %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Biel González, Júlia Soliva  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% AMVO - Course 2025-2026 %%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------%

% -------- PART 2: PRANDTL’S LIFTING LINE MODEL ------------------ %
% -------- APPLIED TO COMPOUND WINGS OF LARGE ASPECT RATIO ------- %

clear; close all; clc;
%% Dades
b  = 15.0;
cr = 0.95;
ct = 0.55;
lambda = ct/cr;
S = (cr+ct)/2*b;
c_bar = (2/3)*cr*(1 + lambda + lambda^2)/(1 + lambda); % MAC [m]
AR = b^2/S;

% Part 1 per HQ300
Cla = 6.935845305522705;   
Cl0 = 0.623151006958811;  

Cd_ala = @(Cl) 0.0183*Cl.^2 - 0.0302*Cl + 0.0187;
alpha_root_deg = 4.0;
alpha_root     = deg2rad(4);

%% Discretització de l'ala
N    = 100; % Per banda
Ntot = 2*N;

y_edges = linspace(-b/2, b/2, Ntot+1);
y_mid   = (y_edges(1:end-1) + y_edges(2:end)) / 2;
dy      = y_edges(2) - y_edges(1);

% Corda local (variació lineal de cr (y=0) a ct (y=+/-b/2))
c_mid = cr + (ct - cr) * (2*abs(y_mid)/b);

% x del quart de corda (ala recta, sense fletxa > x=0 constant)
x_mid = zeros(1, Ntot);

%% Matriu d'influència
Vij = zeros(Ntot, Ntot);
for i = 1:Ntot
    yP = y_mid(i);
    for j = 1:Ntot
        yL = y_edges(j);      % banda esquerra
        yR = y_edges(j+1);    % banda dreta
        dR = yP - yR;
        dL = yP - yL;
        wR = (abs(dR) > 1e-14) *  1/(4*pi*dR);
        wL = (abs(dL) > 1e-14) * -1/(4*pi*dL);

        Vij(i,j) = wR + wL;
    end
end

%%  Valors de twist
twist_values_deg = [0, -1, -2, -3, -4, -5, -6, -7, -8];
n_twists = length(twist_values_deg);

CL_dist    = zeros(1, n_twists); CDind_dist = zeros(1, n_twists); CDv_dist   = zeros(1, n_twists);
Cl_total  = zeros(Ntot, n_twists); Cdi_total = zeros(Ntot, n_twists); Cdv_total = zeros(Ntot, n_twists);
alpha_ind_total  = zeros(Ntot, n_twists); G_total   = zeros(Ntot, n_twists); Cd_total = zeros(Ntot, n_twists);

Qinfmod = 1;  

for jj = 1:n_twists
    theta_tip = deg2rad(twist_values_deg(jj));
    theta_mid = theta_tip * (2*abs(y_mid)/b);

    % Segons pag 28 power. 
    half_c_Cla = 0.5 * c_mid * Cla / Qinfmod;  
    A_mat = eye(Ntot) - bsxfun(@times, half_c_Cla', Vij);
    b_vec = (0.5 * Qinfmod * c_mid .* (Cl0 + Cla*(alpha_root + theta_mid)))';

    Gamma = A_mat \ b_vec;

    % Angle d'atac induït. No afegim k_inf Vij només té en compte la component w. 
    w_ind   = Vij * Gamma;                    
    alpha_ind_loc  = w_ind / Qinfmod;

    % Locals
    Cl_loc  = 2*Gamma ./ (Qinfmod * c_mid');                     
    Cdi_loc = -Cl_loc .* alpha_ind_loc;            
    Cdv_loc = Cd_ala(Cl_loc);

    % Globals
    CL_val    =  2*sum(Gamma * dy) / (Qinfmod * S);
    CDind_val = -2*sum(Gamma .* alpha_ind_loc * dy) / (Qinfmod * S);
    CDv_val   =  sum(Cdv_loc .* c_mid' * dy) / S;
    CD_val = CDind_val + CDv_val;

    CL_dist(jj) = CL_val;
    CDind_dist(jj) = CDind_val;
    CDv_dist(jj) = CDv_val;
    Cl_total(:,jj) = Cl_loc;
    Cdi_total(:,jj) = Cdi_loc;
    Cdv_total(:,jj) = Cdv_loc;
    alpha_ind_total(:,jj)  = alpha_ind_loc;
    G_total(:,jj) = Gamma;

    fprintf('theta_t = %+.0f° | CL = %.4f | CDind = %.5f | CDv = %.5f | CD = %.5f | L/D = %.2f\n', ...
        twist_values_deg(jj), CL_val, CDind_val, CDv_val, CDind_val+CDv_val, ...
        CL_val/(CDind_val+CDv_val));
end

%% Twist amb mínim CDind
[~, idx_opt] = min(CDind_dist);
theta_cd_min_deg = twist_values_deg(idx_opt);
fprintf('\n>>> Twist que minimitza CDind: theta_t = %+.0f°\n\n', theta_cd_min_deg);
eta = y_mid / (b/2);   % coordenada normalitzada

%% Taula
fprintf('\n=======================================================================\n');
fprintf('  theta_t [°]   CL        CDind      CDv        CDtotal     L/D\n');
fprintf('-----------------------------------------------------------------------\n');
theta_max_L_D = 0; valor_max = 0;
for jj = 1:n_twists
    CD = CDind_dist(jj) + CDv_dist(jj);
    fprintf('  %+6.1f       %.5f   %.6f   %.6f   %.6f   %.3f\n', twist_values_deg(jj), CL_dist(jj), CDind_dist(jj), CDv_dist(jj), CD, CL_dist(jj)/CD);
    valor = CL_dist(jj)/CD;
    if valor > valor_max
        valor_max = valor;
        theta_max_L_D = twist_values_deg(jj);
    end
end
fprintf('=======================================================================\n');
fprintf('Twist òptim (maximització L/D): theta_t = %+.0f°\n', theta_max_L_D);

%% Figures
figure;
cmap = parula(n_twists);
hold on
for jj = 1:n_twists
    plot(eta, Cl_total(:,jj), 'Color', cmap(jj,:), 'LineWidth', 0.5, 'DisplayName', sprintf('\\theta_t = %+.0f°', twist_values_deg(jj)));
end
xlabel('2y/b','FontSize',12)
ylabel('C_l','FontSize',12)
title('Spanwise distribution of section lift coefficient','FontSize',12)
legend('Location','south','NumColumns',3,'FontSize',9)
grid on; xlim([-1 1])

figure;
cmap = parula(n_twists);
hold on
for jj = 1:n_twists
    plot(eta, Cdi_total(:,jj), 'Color', cmap(jj,:), 'LineWidth', 0.5, 'DisplayName', sprintf('\\theta_t = %+.0f°', twist_values_deg(jj)));
end
xlabel('2y/b','FontSize',12)
ylabel('C_di','FontSize',12)
title('Spanwise distribution of section induced drag coefficient','FontSize',12)
legend('Location','south','NumColumns',3,'FontSize',9)
grid on; xlim([-1 1])

figure;
subplot(1,2,1)
plot(twist_values_deg, CL_dist, 'bo-','LineWidth',1,'MarkerFaceColor','b','MarkerSize',3)
xlabel('\theta_t [°]','FontSize',11)
ylabel('C_L','FontSize',11)
title('Total Lift','FontSize',11)
grid on; grid minor
subplot(1,2,2)
plot(twist_values_deg, CDind_dist, 'r^-','LineWidth',1,'MarkerFaceColor','r','MarkerSize',3,'DisplayName','C_{D,ind}')
hold on
plot(twist_values_deg, CDind_dist+CDv_dist, 'ks-','LineWidth',1,'MarkerFaceColor','k','MarkerSize',3,'DisplayName','C_{D,total}')
xline(theta_cd_min_deg,'--','Color',[0 0.6 0],'LineWidth',1.5,'Label',sprintf('\\theta_t=%+.0f°',theta_cd_min_deg))
xlabel('\theta_t [°]','FontSize',11)
ylabel('C_D','FontSize',11)
title('Total Drag','FontSize',11)
legend('Location','best','FontSize',9); grid on; grid minor
sgtitle(sprintf('Wing twist effect | \\alpha = %.0f°', alpha_root_deg), 'FontSize',12,'FontWeight','bold')