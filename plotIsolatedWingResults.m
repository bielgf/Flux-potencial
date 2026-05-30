function plotIsolatedWingResults(twist_val,eta,Cl_vec,Cd_ind_vec_1,Cd_visc_vec_1,alpha_ind_vec,CL,CD_ind,CD,theta_cd_min_deg,theta_max_L_D,alpha)

twist_val_filt     = twist_val(mod(rad2deg(twist_val),1) == 0);
Cl_vec_filt        = Cl_vec(:,mod(rad2deg(twist_val),1) == 0);
Cd_ind_vec_1_filt  = Cd_ind_vec_1(:,mod(rad2deg(twist_val),1) == 0);
Cd_visc_vec_1_filt = Cd_visc_vec_1(:,mod(rad2deg(twist_val),1) == 0);
alpha_ind_vec_filt = alpha_ind_vec(:,mod(rad2deg(twist_val),1) == 0);
CL_filt            = CL(mod(rad2deg(twist_val),1) == 0);
CD_ind_filt        = CD_ind(mod(rad2deg(twist_val),1) == 0);
CD_filt            = CD(mod(rad2deg(twist_val),1) == 0);

n_twist_filt = length(twist_val_filt);

% ── Spanwise lift ────
figure;
cmap = parula(n_twist_filt);
hold on
for jj = 1:n_twist_filt
    plot(eta, Cl_vec_filt(:,jj), 'Color', cmap(jj,:), 'LineWidth', 0.5, 'DisplayName', sprintf('\\theta_t = %+.0f°', rad2deg(twist_val_filt(jj))));
end
xlabel('2y/b','FontSize',12)
ylabel('C_l','FontSize',12)
title('Spanwise distribution of section lift coefficient','FontSize',12)
legend('Location','south','NumColumns',3,'FontSize',9)
grid on; xlim([-1 1])

% ── Spanwise total drag ───
figure;
cmap = parula(n_twist_filt);
hold on
for jj = 1:n_twist_filt
    plot(eta, Cd_ind_vec_1_filt(:,jj), 'Color', cmap(jj,:), 'LineWidth', 0.5, 'DisplayName', sprintf('\\theta_t = %+.0f°', rad2deg(twist_val_filt(jj))));
end
xlabel('2y/b','FontSize',12)
ylabel('C_di','FontSize',12)
title('Spanwise distribution of section induced drag coefficient','FontSize',12)
legend('Location','south','NumColumns',3,'FontSize',9)
grid on; xlim([-1 1])

% ── Spanwise viscous drag ────
figure;
cmap = parula(n_twist_filt);
hold on
for jj = 1:n_twist_filt
    plot(eta, Cd_visc_vec_1_filt(:,jj), 'Color', cmap(jj,:), 'LineWidth', 0.5, ...
        'DisplayName', sprintf('\\theta_t = %+.0f°', rad2deg(twist_val_filt(jj))));
end
xlabel('2y/b', 'FontSize', 12)
ylabel('C_{d,visc}', 'FontSize', 12)
title('Spanwise distribution of section viscous drag coefficient', 'FontSize', 12)
legend('Location', 'south', 'NumColumns', 3, 'FontSize', 9)
grid on; xlim([-1 1])

% ── Spanwise induced angle of attack ───
figure;
cmap = parula(n_twist_filt);
hold on
for jj = 1:n_twist_filt
    plot(eta, rad2deg(alpha_ind_vec_filt(:,jj)), 'Color', cmap(jj,:), 'LineWidth', 0.5, ...
        'DisplayName', sprintf('\\theta_t = %+.0f°', rad2deg(twist_val_filt(jj))));
end
xlabel('2y/b', 'FontSize', 12)
ylabel('\alpha_{ind} [°]', 'FontSize', 12)
title('Spanwise distribution of induced angle of attack', 'FontSize', 12)
legend('Location', 'south', 'NumColumns', 3, 'FontSize', 9)
grid on; xlim([-1 1])

figure;
subplot(1,2,1)
plot(rad2deg(twist_val_filt), CL_filt, 'bo-','LineWidth',1,'MarkerFaceColor','b','MarkerSize',3)
xlabel('\theta_t [°]','FontSize',11)
ylabel('C_L','FontSize',11)
title('Total Lift','FontSize',11)
grid on; grid minor
subplot(1,2,2)
plot(rad2deg(twist_val_filt), CD_ind_filt, 'r^-','LineWidth',1,'MarkerFaceColor','r','MarkerSize',3,'DisplayName','C_{D,ind}')
hold on
plot(rad2deg(twist_val_filt), CD_filt, 'ks-','LineWidth',1,'MarkerFaceColor','k','MarkerSize',3,'DisplayName','C_{D,total}')
xline(theta_cd_min_deg,'--','Color',[0 0.6 0],'LineWidth',1.5,'DisplayName','Minimum C_{d_{ind}}','Label',sprintf('\\theta_t=%+.2f°',theta_cd_min_deg))
xline(theta_max_L_D,'--','Color',[0 0 0.6],'LineWidth',1.5,'DisplayName','Maximum C_L/C_D','Label',sprintf('\\theta_t=%+.2f°',theta_max_L_D))
xlabel('\theta_t [°]','FontSize',11)
ylabel('C_D','FontSize',11)
title('Total Drag','FontSize',11)
legend('Location','best','FontSize',9); grid on; grid minor
sgtitle(sprintf('Wing twist effect | \\alpha = %.0f°', rad2deg(alpha)), 'FontSize',12,'FontWeight','bold')

figure;
plot(rad2deg(twist_val), CL./CD, 'bo-', 'LineWidth', 1,'MarkerFaceColor', 'b', 'MarkerSize', 2)
xline(theta_max_L_D, '--', 'Color', 'k', 'LineWidth', 1.5, 'Label', sprintf('\\theta_t = %.2f°', theta_max_L_D), 'LabelVerticalAlignment', 'middle');
xlabel('\theta_t [°]', 'FontSize', 12)
ylabel('C_L / C_D', 'FontSize', 12)
title('Lift-to-drag ratio vs wing tip twist | \alpha = 4°', 'FontSize', 12)
grid on; grid minor


end