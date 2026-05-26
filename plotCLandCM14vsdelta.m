function plotCLandCM14vsdelta(delta_e,CL_table,CM14_table)

    figure('Name','Cl Slope','Color','w');
    hold on; grid on;
    plot(rad2deg(delta_e),CL_table(:,end),'LineWidth',1.5);
    title('LIFT COEFFICIENT SLOPE (C_l vs \delta_e for \alpha = 4{\circ})','FontName','Times New Roman','FontWeight','normal');
    xlabel(sprintf('\\delta_e (degrees)'),'FontSize', 12);
    ylabel('C_l','FontSize', 12);

    figure('Name','Cm14 Slope','Color','w');
    hold on; grid on;
    plot(rad2deg(delta_e),CM14_table(:,end),'LineWidth',1.5);
    title('PITCHING MOMENT COEFFICIENT SLOPE (C_{m_{1/4}} vs \delta_e for \alpha = 4{\circ})','FontName','Times New Roman','FontWeight','normal');
    xlabel(sprintf('\\delta_e (degrees)'),'FontSize', 12);
    ylabel('C_{m_{1/4}}','FontSize', 12);

end