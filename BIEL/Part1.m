clear;clc;close all

%% Preprocess

% Read the data file
Ndiv = 16;
filename = sprintf('C:\\Users\\Biel\\Desktop\\UNI\\MUEA\\Q2\\AMVO\\POTENCIAL\\Flux potencial\\HQ_300\\HQ300_%.0f.txt',Ndiv);

% Open and read the file
data = load(filename);  % or use importdata/readmatrix

% Extract variables
index = data(:, 1);
X     = [data(:, 2) data(:, 3)];

c      = 1;

Nc = zeros(Ndiv,2);
Tc = zeros(Ndiv,2);

ca = zeros(Ndiv,1);
sa = zeros(Ndiv,1);
l  = zeros(Ndiv,1);

Xc    = zeros(Ndiv,2);
delta = zeros(Ndiv,2);

cl    = zeros(Ndiv,1);
cp    = zeros(Ndiv,1);
cm1_4 = zeros(Ndiv,1);

for jj = 1:Ndiv
    l(jj)       = sqrt((X(jj,1) - X(jj+1,1))^2 + (X(jj,2) - X(jj+1,2))^2);
    Xc(jj,:)    = (X(jj,:) + X(jj+1,:))/2;
    delta(jj,:) = X(jj+1,:) - X(jj,:);
    ca(jj)      = (X(jj+1,1) - X(jj,1))/l(jj);
    sa(jj)      = (X(jj,2) - X(jj+1,2))/l(jj);

    Nc(jj,:)  = [sa(jj,1),ca(jj,1)];
    Tc(jj,:)  = [ca(jj,1),-sa(jj,1)];
end

%% Process

a = zeros(Ndiv,Ndiv);
b = zeros(Ndiv,1);

alpha = deg2rad(0);

rho = 1.225;

Qinfmod = 1;
Qinf    = Qinfmod*[cos(alpha),sin(alpha)];

for ii = 1:Ndiv
    b(ii) = -dot(Qinf,Tc(ii,:));
    for jj = 1:Ndiv
        if jj ~= ii
            Xc_ij_pan = (Xc(ii,1) - X(jj,1))*ca(jj) - (Xc(ii,2) - X(jj,2))*sa(jj);
            Zc_ij_pan = (Xc(ii,1) - X(jj,1))*sa(jj) + (Xc(ii,2) - X(jj,2))*ca(jj);

            r1 = sqrt(Xc_ij_pan^2 + Zc_ij_pan^2);
            r2 = sqrt((Xc_ij_pan - l(jj))^2 + Zc_ij_pan^2);

            theta1 = atan2(Zc_ij_pan,Xc_ij_pan);
            theta2 = atan2(Zc_ij_pan,Xc_ij_pan - l(jj));

            u_ij_pan = (theta2 - theta1)/(2*pi);
            w_ij_pan = (1/(4*pi))*log(r2^2/r1^2);

            u_ij = u_ij_pan*ca(jj) + w_ij_pan*sa(jj);
            w_ij = -u_ij_pan*sa(jj) + w_ij_pan*ca(jj);

            a(ii,jj) = dot([u_ij,w_ij],Tc(ii,:));
        else
            a(ii,jj) = -0.5;
        end
    end
end

K         = fix(Ndiv/4);
a(K,:)    = 0;
a(K,1)    = 1;
a(K,Ndiv) = 1;
b(K)      = 0;

gamma = a\b;
gamma(K) = 0.5*(gamma(K-1)+gamma(K+1));

for ii=1:Ndiv
     cl(ii)    = 2*gamma(ii)*l(ii)/Qinfmod;
     cp(ii)    = 1 - (gamma(ii)/Qinfmod)^2;
     cm1_4(ii) = cp(ii)*((Xc(ii,1)/c)*(delta(ii,1)/c) + (Xc(ii,2)/c)*(delta(ii,2)/c)) - 0.25*cl(ii);
end

CL    = sum(cl);
L     = CL*0.5*Qinfmod^2*2*c*rho;
CM1_4 = sum(cm1_4);
M1_4  = CM1_4*0.5*rho*Qinfmod^2*c^2;

%% Kármán-Tsien Compressibility Correction & Critical Mach Number

gam = 1.4;

Cp0     = cp;
Cp0_min = min(Cp0);

% Kármán-Tsien correction
KT = @(Cp0_val, M) Cp0_val./(sqrt(1 - M.^2) + (M.^2./(1 + sqrt(1 - M.^2))).*(Cp0_val/2));

% Critical Cp (local Mach = 1, isentropic)
Cp_crit = @(M) (2./(gam.*M.^2)).*(((2/(gam+1)).*(1 + (gam-1)/2.*M.^2)).^(gam/(gam-1)) - 1 );

residual = @(M) KT(Cp0_min, M) - Cp_crit(M);

% Bracket search
% KT is only physically valid where its denominator > 0
% For very negative Cp0, the denominator flips sign at high M → non-physical
M_scan   = linspace(0.02, 0.99, 5000);
denom_KT = sqrt(1 - M_scan.^2) + (M_scan.^2./(1 + sqrt(1 - M_scan.^2))).*(Cp0_min/2);

valid    = denom_KT > 0;                  % physical region only
M_valid  = M_scan(valid);
res_scan = arrayfun(residual,M_valid);

% Find first sign change
sc = find(diff(sign(res_scan)) ~= 0, 1, 'first');

if isempty(sc)
    error('No critical Mach found. Cp0_min = %.4f — body may have no suction.', Cp0_min);
end

Mcr = fzero(residual, [M_valid(sc), M_valid(sc+1)]);
fprintf('=== Critical Mach Number: Mcr = %.4f ===\n', Mcr);

%% Post-processing Visualization

theta_plot = ((1:Ndiv)' - 0.5)*360/Ndiv;

% Analytic solutions (clockwise convention to match panel ordering)
theta_an_rad = deg2rad(linspace(0,360,1000));
Vt_an        = -2*Qinfmod*sin(theta_an_rad + alpha);            % tangential velocity
Cp_an        =  1 - (2*sin(theta_an_rad + alpha)).^2;           % pressure coefficient

% ── Figure 1: Vortex strength distribution ─────────────────────────────────
figure('Name','Vortex Strength Distribution','Color','w');
hold on; axis equal off;

plot(X(:,1), X(:,2), 'k-', 'LineWidth', 1.5);

% sc1 = 0.25 / max(abs(gamma));               % scale so max height = 0.25
sc1 = 0.06;
for jj = 1:Ndiv
    col = 'b';
    gammaPlot = gamma(jj);
    if gamma(jj) < 0
        col = 'r';
        gammaPlot = -gamma(jj);
    end
    p1 = X(jj,:);
    p2 = X(jj+1,:);
    p3 = p2 + gammaPlot*sc1*Nc(jj,:);      % outer corners
    p4 = p1 + gammaPlot*sc1*Nc(jj,:);
    fill([p1(1),p2(1),p3(1),p4(1)], ...
         [p1(2),p2(2),p3(2),p4(2)], col, ...
         'EdgeColor',col, 'LineWidth',0.3, 'FaceAlpha',0.5);
end

text(0.50,  0.25, '\gamma > 0', 'Color','b', 'FontSize',12, ...
      'FontName','Times New Roman');
text(0.50, -0.15, '\gamma < 0', 'Color','r', 'FontSize',12, ...
      'FontName','Times New Roman');
text(0.50,  0.04, sprintf('N_{div} = %d', Ndiv), 'FontSize',11, ...
      'FontName','Times New Roman', 'HorizontalAlignment','center');
title('DISTRIBUTION OF VORTEX STRENGTH', ...
      'FontName','Times New Roman', 'FontWeight','normal');

% ── Figure 2: Tangential velocity distribution ─────────────────────────────
figure('Name','Surface Velocity Distribution','Color','w');
hold on; grid on; box on; grid minor;

plot(theta_plot,  gamma,  'b-', 'LineWidth',1.5, 'DisplayName','CSV');
plot(rad2deg(theta_an_rad), Vt_an, 'bo', 'MarkerSize',4,  'DisplayName','Analytic');

xlim([0 360]); xticks(0:50:350);
xlabel('\theta (deg.)',  'FontName','Times New Roman', 'FontSize',12);
ylabel('V_t',           'FontName','Times New Roman', 'FontSize',12);
title({'DISTRIBUTION OF VELOCITY ON THE SURFACE','OF THE CYLINDER'}, ...
       'FontName','Times New Roman', 'FontWeight','normal');
legend('Location','southeast', 'FontSize',10);

ysc = max(abs(gamma));
text(130,  0.55*ysc, sprintf('\\alpha = %g deg.', rad2deg(alpha)), ...
     'FontSize',11, 'FontName','Times New Roman');
text(130,  0.25*ysc, sprintf('N_{div} = %d', Ndiv), ...
     'FontSize',11, 'FontName','Times New Roman');

% ── Figure 3: Cp line plot ──────────────────────────────────────────────────
figure('Name','Cp Distribution','Color','w');
hold on; grid on; box on; grid minor;

plot(X(:,1), X(:,2), 'k-', 'LineWidth', 1.5);

xlim([0 360]); xticks(0:50:350);
xlabel('\theta (deg.)', 'FontName','Times New Roman', 'FontSize',12);
ylabel('C_p',           'FontName','Times New Roman', 'FontSize',12);
title('DISTRIBUTION OF PRESSURE COEFFICIENT', ...
      'FontName','Times New Roman', 'FontWeight','normal');

cp_bot = min(cp);                    % most negative Cp, for text placement
text(130, 0.30*cp_bot, sprintf('\\alpha = %g deg.', rad2deg(alpha)), ...
     'FontSize',11, 'FontName','Times New Roman');
text(130, 0.55*cp_bot, sprintf('N_{div} = %d', Ndiv), ...
     'FontSize',11, 'FontName','Times New Roman');

% ── Figure 4: Cp vectors on cylinder body (with arrowheads) ────────────────
figure('Name','Cp on Body','Color','w');
hold on; axis equal off;

plot(X(:,1), X(:,2), 'k-', 'LineWidth', 1.5);

sc4 = 0.06 / max(abs(cp));   % scale: max vector length = 0.3
hs  = 0.4;                  % MaxHeadSize (relative to arrow length)

mask_suc  = cp <  0;        % suction panels
mask_pres = cp >= 0;        % pressure panels

% ── Suction (Cp < 0, red): tail AT surface, head OUTSIDE, points outward ──
x0_r = Xc(mask_suc,1);
y0_r = Xc(mask_suc,2);
dx_r = abs(cp(mask_suc)).*sc4.*Nc(mask_suc,1);
dy_r = abs(cp(mask_suc)).*sc4.*Nc(mask_suc,2);
quiver(x0_r, y0_r, dx_r, dy_r, 0, 'r', ...
       'MaxHeadSize', hs, 'LineWidth', 0.8);

% ── Pressure (Cp > 0, blue): tail OUTSIDE, head AT surface, points inward ──
x0_b = Xc(mask_pres,1) + cp(mask_pres).*sc4.*Nc(mask_pres,1);
y0_b = Xc(mask_pres,2) + cp(mask_pres).*sc4.*Nc(mask_pres,2);
dx_b = Xc(mask_pres,1) - x0_b;
dy_b = Xc(mask_pres,2) - y0_b;
quiver(x0_b, y0_b, dx_b, dy_b, 0, 'b', ...
       'MaxHeadSize', hs, 'LineWidth', 0.8);

% ── Annotations ────────────────────────────────────────────────────────────
x_txt = R + max(abs(cp))*sc4 + 0.12;
text(x_txt,  0.15, sprintf('\\alpha = %g deg.', rad2deg(alpha)), ...
     'FontSize',11, 'FontName','Times New Roman');
text(x_txt,  0.00, sprintf('C_l = %.4f', CL), ...
     'FontSize',11, 'FontName','Times New Roman');
text(x_txt, -0.15, sprintf('N_{div} = %d', Ndiv), ...
     'FontSize',11, 'FontName','Times New Roman');
title('DISTRIBUTION OF PRESSURE COEFFICIENT', ...
      'FontName','Times New Roman', 'FontWeight','normal');

% Pel segon exercici part del treball - Part 1
% A = [A11 A12 ; A21 A22]

% Hi han dues condicions de Kutta, ens hem de carregar dues files de la
% matriu A (corresponent a panells del intradós dels perfils corresponents)