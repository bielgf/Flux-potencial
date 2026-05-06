function plotGammaDistribution(Ndiv,X,Nc,gamma)

    figure('Name','Vortex Strength Distribution','Color','w');
    hold on; axis equal off;
    
    plot(X(:,1), X(:,2), 'k-', 'LineWidth', 1.5);
    
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
        p3 = p2 + gammaPlot*sc1*Nc(jj,:);
        p4 = p1 + gammaPlot*sc1*Nc(jj,:);
        fill([p1(1),p2(1),p3(1),p4(1)],[p1(2),p2(2),p3(2),p4(2)], col,'EdgeColor',col,'LineWidth',0.3,'FaceAlpha',0.5);
    end
    
    text(0.50,0.25,'\gamma > 0','Color','b','FontSize',12,'FontName','Times New Roman');
    text(0.50,-0.15,'\gamma < 0','Color','r','FontSize',12,'FontName','Times New Roman');
    text(0.50,0.04,sprintf('N_{div} = %d',Ndiv),'FontSize',11,'FontName','Times New Roman','HorizontalAlignment','center');
    title('DISTRIBUTION OF VORTEX STRENGTH','FontName','Times New Roman','FontWeight','normal');

end