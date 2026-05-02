% =========================================================================
% AERODINÀMICA, MECÀNICA DE VOL I ORBITAL (2025-2026)
% EXERCICI 2: Flux potencial 2D al voltant d'un cilindre (CSV)
% =========================================================================

clear; clc; close all;

% VARIABLES

Rc   = 0.5;       % Radi del cilindre (corda c = 2*Rc = 1)
Qinf = 1.0;
rho  = 1.225;

%% 1. DISTRIBUCIÓ D'INTENSITAT DE VÒRTEX

N = 32;
alpha_deg = 0;
gamma = solve_csv_cylinder(N, Rc, Qinf, alpha_deg, rho);

% Geometria
theta_nodes = linspace(2*pi, 0, N+1)';
X = Rc * cos(theta_nodes);
Z = Rc * sin(theta_nodes);
xc = (X(1:end-1) + X(2:end)) / 2;
zc = (Z(1:end-1) + Z(2:end)) / 2;

% Figura 1: Distribució de Vòrtexs
figure('Name', 'Vortex Strength', 'Color', 'w');
hold on; axis equal; grid on;
title(sprintf('DISTRIBUTION OF VORTEX STRENGTH\nN_{div} = %d', N));
plot(X, Z, 'k', 'LineWidth', 1);

% Rectangles d'intensitat
scale = 0.15; 
for i = 1:N
    % Vector normal 
    dx = X(i+1) - X(i);
    dz = Z(i) - Z(i+1);
    l = sqrt(dx^2 + dz^2);
    nx = dz/l; 
    nz = dx/l;
    
    h = abs(gamma(i) * scale);
    color_g = 'b'; if gamma(i) < 0; color_g = 'r'; end
    
    p1 = [X(i), Z(i)];
    p2 = [X(i+1), Z(i+1)];
    p3 = p2 + [nx, nz]*h;
    p4 = p1 + [nx, nz]*h;
    
    fill([p1(1) p2(1) p3(1) p4(1)], [p1(2) p2(2) p3(2) p4(2)], color_g, 'FaceAlpha', 0.1, 'EdgeColor', color_g);
end
text(0.6, 0.4, '\gamma > 0', 'Color', 'b', 'FontSize', 12);
text(0.6, -0.4, '\gamma < 0', 'Color', 'r', 'FontSize', 12);

%% 2. DISTRIBUCIÓ DE VELOCITAT I CP 

N = 256;
alpha_deg = 0;
[gamma, xc, zc, Vt, Cp] = solve_csv_cylinder(N, Rc, Qinf, alpha_deg, rho);

theta_c = atan2(zc, xc);
theta_c_plot = mod(-theta_c * 180/pi, 360); % negative bc clockwise 
[theta_c_plot, sort_idx] = sort(theta_c_plot);
Cp = Cp(sort_idx);
Vt = Vt(sort_idx);

alpha_rad = alpha_deg * pi / 180;
Vt_ana = -2*Qinf*sin(theta_c(sort_idx) - alpha_rad) - 2*Qinf*sin(alpha_rad);
Cp_ana = 1 - (Vt_ana / Qinf).^2;

figure('Name', 'Velocity and Cp', 'Color', 'w', 'Position', [100 100 600 800]);

subplot(2,1,1); hold on; grid on;
title('DISTRIBUTION OF VELOCITY ON THE SURFACE OF THE CYLINDER');
plot(theta_c_plot, Vt, 'c-', 'LineWidth', 1, 'DisplayName', 'CSV');
plot(theta_c_plot, Vt_ana, 'bo', 'MarkerSize', 5, 'DisplayName', 'Analytic');
xlabel('\theta (deg.)'); ylabel('V_t');
legend('Location', 'best'); xlim([0 360]);

subplot(2,1,2); hold on; grid on;
title('DISTRIBUTION OF PRESSURE COEFFICIENT');
plot(theta_c_plot, Cp, 'c-', 'LineWidth', 1, 'DisplayName', 'CSV');
plot(theta_c_plot, Cp_ana, 'b-', 'LineWidth', 1, 'DisplayName', 'Analytic');
xlabel('\theta (deg.)'); ylabel('C_p');
xlim([0 360]); ylim([-3 1]);

%% 3. CONVERGENCIA DEL COEFICIENT DE SUSTENTACIÓ 
alpha_deg = 6;
alpha_rad = alpha_deg * pi / 180;
N_list = [32, 64, 128, 256, 512, 1024];
Cl_csv = zeros(size(N_list));
Cl_ana = 4 * pi * sin(alpha_rad); % Fòrmula analítica

for k = 1:length(N_list)
    [~, ~, ~, ~, ~, Cl] = solve_csv_cylinder(N_list(k), Rc, Qinf, alpha_deg, rho);
    Cl_csv(k) = Cl;
end

figure('Name', 'Cl Convergence', 'Color', 'w');
hold on; grid on;
title('LIFT COEFFICIENT CONVERGENCE');
plot(N_list, Cl_csv, 'bo', 'MarkerSize', 8, 'DisplayName', 'CSV');
yline(Cl_ana, 'k--', 'DisplayName', 'Analytic');
xlabel('N_{div}'); ylabel('C_l');
legend('Location', 'southeast');
text(700, 1.28, '\alpha = 6 deg.', 'FontSize', 12);

%% 4. DISTRIBUCIÓ DE CP PER DIFERENTS ANGLES 
alphas_test = [0, 2, 4, 6, 8];
N = 256;

figure('Name', 'Cp Multiple Alphas', 'Color', 'w');
hold on; grid on;
title('DISTRIBUTION OF PRESSURE COEFFICIENT');
colors = lines(length(alphas_test));

for k = 1:length(alphas_test)
    [~, xc, zc, ~, Cp] = solve_csv_cylinder(N, Rc, Qinf, alphas_test(k), rho);
    theta_c = atan2(zc, xc);
    theta_c_plot = mod(-theta_c * 180/pi, 360); % negative bc clockwise 
    [theta_c_plot, sort_idx] = sort(theta_c_plot);
    Cp = Cp(sort_idx);
    
    plot(theta_c_plot, Cp, 'Color', colors(k,:), 'LineWidth', 1, ...
        'DisplayName', sprintf('\\alpha = %d deg.', alphas_test(k)));
end
xlabel('\theta (deg.)'); ylabel('C_p');
legend('Location', 'southwest'); xlim([0 360]);

%% 5. COEFICIENT DE SUSTENTACIÓ VS ANGLE D'ATAC 
alpha_sweep = 0:1:8;
Cl_sweep_csv = zeros(size(alpha_sweep));
Cl_sweep_ana = 4 * pi * sin(alpha_sweep * pi / 180);

for k = 1:length(alpha_sweep)
    [~, ~, ~, ~, ~, Cl] = solve_csv_cylinder(N, Rc, Qinf, alpha_sweep(k), rho);
    Cl_sweep_csv(k) = Cl;
end

figure('Name', 'Cl vs Alpha', 'Color', 'w');
hold on; grid on;
title('LIFT COEFFICIENT VS. ANGLE OF ATTACK');
plot(alpha_sweep, Cl_sweep_csv, 'k-', 'LineWidth', 1, 'DisplayName', 'CSV');
plot(alpha_sweep, Cl_sweep_ana, 'r:', 'LineWidth', 1.5, 'DisplayName', 'Analytic');
xlabel('\alpha (deg.)'); ylabel('C_l');
legend('Location', 'southeast');


%% ========================================================================
% SOLVER MÈTODE DE PANELLS (CSV)
% =========================================================================
function [gamma, xc, zc, Vt, Cp, Cl] = solve_csv_cylinder(N, Rc, Qinf, alpha_deg, rho)
    
    % 1. Discretizació de la superficie (Clockwise)
    theta_nodes = linspace(2*pi, 0, N+1)';
    X = Rc * cos(theta_nodes);
    Z = Rc * sin(theta_nodes);
    
    % 2. Initzialització de vectors
    xc = zeros(N,1); zc = zeros(N,1);
    l = zeros(N,1); C = zeros(N,1); S = zeros(N,1);
    
    for i = 1:N
        xc(i) = (X(i) + X(i+1))/2;
        zc(i) = (Z(i) + Z(i+1))/2;
        l(i) = sqrt((X(i+1)-X(i))^2 + (Z(i)-Z(i+1))^2);
        
        C(i) = (X(i+1)-X(i))/l(i); 
        S(i) = (Z(i)-Z(i+1))/l(i); 
    end
    
    % 3. Matriu de Coeficients d'Influència (A) i Terme Independent (b)
    alpha_rad = alpha_deg * pi / 180;
    Qx = Qinf * cos(alpha_rad);
    Qz = Qinf * sin(alpha_rad);
    
    A = zeros(N,N);
    b = zeros(N,1);
    
    for i = 1:N
        b(i) = -(Qx * C(i) - Qz * S(i)); 
        
        for j = 1:N
            if i == j
                A(i,j) = -0.5; % Auto-inducció
            else
                % Coordenades del punt de control i en el sistema local del panell j
                dx = xc(i) - X(j);
                dz = zc(i) - Z(j);
                x_local = dx * C(j) - dz * S(j);
                z_local = dx * S(j) + dz * C(j);
                
                % Radis y angles
                r1_sq = x_local^2 + z_local^2;
                r2_sq = (x_local - l(j))^2 + z_local^2;
                theta1 = atan2(z_local, x_local);
                theta2 = atan2(z_local, x_local - l(j));
                
                % Velocitats induïdes
                u_local = (theta2 - theta1) / (2*pi);
                w_local = (1/(4*pi)) * log(r2_sq / r1_sq);
                
                % Canvi a sistema global
                u_glob = u_local * C(j) + w_local * S(j);
                w_glob = -u_local * S(j) + w_local * C(j);
                
                % Producte escalar amb vector tangent
                A(i,j) = u_glob * C(i) - w_glob * S(i);
            end
        end
    end
    
    % 4. Condició de Kutta
    iel = round(N/4); 
    A(iel, :) = 0;
    A(iel, 1) = 1;
    A(iel, N) = 1;
    b(iel) = 0;
    
    % 5. Resoldre sistema
    gamma = A \ b;
    
    if iel > 1 && iel < N
        gamma(iel) = (gamma(iel-1) + gamma(iel+1))/2;
    end
    
    % 6. Post-Processing
    Vt = -gamma; % Velocidad tangencial
    Cp = 1 - (Vt / Qinf).^2;
    
    % Coeficient de sustentació (Teorema Kutta-Joukowski)
    L_prime = rho * Qinf * sum(gamma .* l);
    chord = 2*Rc;
    Cl = L_prime / (0.5 * rho * Qinf^2 * chord);
end