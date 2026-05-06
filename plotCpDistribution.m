function plotCpDistribution(X,Nc,Xc,cp,alpha,Ndiv,CL)

    figure('Name','Cp on Body','Color','w');
    hold on; axis equal off;
    
    plot(X(:,1), X(:,2), 'k-', 'LineWidth', 1.5);
    
    sc4 = 0.2/max(abs(cp));      % scale: max vector length = 0.3
    hs  = 0.1;                   % MaxHeadSize (relative to arrow length)
    
    mask_suc  = cp <  0;         % suction panels
    mask_pres = cp >= 0;         % pressure panels
    
    x0_r = Xc(mask_suc,1);
    y0_r = Xc(mask_suc,2);
    dx_r = abs(cp(mask_suc)).*sc4.*Nc(mask_suc,1);
    dy_r = abs(cp(mask_suc)).*sc4.*Nc(mask_suc,2);
    quiver(x0_r,y0_r,dx_r,dy_r,0,'r','MaxHeadSize',hs,'LineWidth',0.8);
    
    x0_b = Xc(mask_pres,1) + cp(mask_pres).*sc4.*Nc(mask_pres,1);
    y0_b = Xc(mask_pres,2) + cp(mask_pres).*sc4.*Nc(mask_pres,2);
    dx_b = Xc(mask_pres,1) - x0_b;
    dy_b = Xc(mask_pres,2) - y0_b;
    quiver(x0_b,y0_b,dx_b,dy_b,0,'b','MaxHeadSize',hs,'LineWidth',0.8);

    text(0.3,0.42, sprintf('\\alpha = %g deg.', rad2deg(alpha)),'FontSize',11, 'FontName','Times New Roman');
    text(0.3,0.35, sprintf('C_l = %.4f', CL),'FontSize',11, 'FontName','Times New Roman');
    text(0.3,0.3, sprintf('N_{div} = %d', Ndiv),'FontSize',11,'FontName','Times New Roman');
    title('DISTRIBUTION OF PRESSURE COEFFICIENT','FontName','Times New Roman','FontWeight','normal');

end