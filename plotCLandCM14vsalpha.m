function plotCLandCM14vsalpha(alpha,CL_table,CM14_table)

    figure('Name','Cl Slope','Color','w');
    hold on; grid on;
    plot(rad2deg(alpha),CL_table(:,end),'LineWidth',1.5);
    title('LIFT COEFFICIENT SLOPE (C_l vs \alpha)','FontName','Times New Roman','FontWeight','normal');
    xlabel(sprintf('\\alpha (degrees)'),'FontSize', 12);
    ylabel('C_l','FontSize', 12);
    Clalpha = polyfit(alpha,CL_table(:,end),1);
    text(4,1,sprintf('C_{l_{\\alpha}} = %.4f',Clalpha(1)))

    figure('Name','Cm14 Slope','Color','w');
    hold on; grid on;
    plot(rad2deg(alpha),CM14_table(:,end),'LineWidth',1.5);
    title('PITCHING MOMENT COEFFICIENT SLOPE (C_{m_{1/4}} vs \alpha)','FontName','Times New Roman','FontWeight','normal');
    xlabel(sprintf('\\alpha (degrees)'),'FontSize', 12);
    ylabel('C_{m_{1/4}}','FontSize', 12);
    Cmalpha = polyfit(alpha,CM14_table(:,end),1);
    text(4,-0.135,sprintf('C_{m_{\\alpha}} = %.4f',Cmalpha(1)))

end