function plotConvergence(Ndiv,CL_table,alpha)

    figure('Name','Cl Convergence','Color','w');
    hold on; grid on;

    for ii = 1:size(CL_table, 1)
        plot(Ndiv, CL_table(ii, :),'-o', 'LineWidth', 1.5, 'DisplayName', sprintf('\\alpha = %g^o', rad2deg(alpha(ii))));
    end
    
    xlabel('N_{div}');
    ylabel('C_l');
    title('LIFT COEFFICIENT CONVERGENCE','FontName','Times New Roman','FontWeight','normal');
    legend('Location', 'best');

end