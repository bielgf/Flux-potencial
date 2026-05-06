function plotCpChordDistribution(X,cp,alpha,Ndiv)

    % Calcular las coordenadas X e Y de los centros de los paneles (puntos de control)
    X_nodes = X(:, 1);
    X_mid   = (X_nodes(1:end-1) + X_nodes(2:end))/2; % Tamaño Ndiv
    
    % Encontrar el índice del Borde de Ataque (Leading Edge)
    % Es el punto donde la coordenada X es mínima
    [~, LE_idx] = min(X_mid);
    
    % Separar los índices en Intradós y Extradós
    % Asumimos el orden estándar de los métodos de paneles: 
    % Comienza en Borde de Salida inferior -> Borde de Ataque -> Borde de Salida superior.
    idx_intrados = 1:LE_idx;
    idx_extrados = LE_idx+1:length(X_mid);

    % Crear la figura
    figure;
    hold on; 
    grid on; % Activar la cuadrícula
    
    % Dibujar la línea negra continua que une todos los puntos del Cp
    plot(X_mid, cp, 'k-', 'LineWidth', 1);
    
    % Dibujar los puntos del Intradós (Círculos azules vacíos)
    p_int = plot(X_mid(idx_intrados), cp(idx_intrados), 'bo', 'MarkerSize', 6);
    
    % Dibujar los puntos del Extradós (Círculos rojos vacíos)
    p_ext = plot(X_mid(idx_extrados), cp(idx_extrados), 'ro', 'MarkerSize', 6);
    
    xlabel('X', 'FontSize', 12);
    ylabel('C_p', 'FontSize', 12);
    legend([p_int, p_ext], {'Intrados', 'Extrados'}, 'Location', 'southeast');
    set(gca, 'YDir', 'reverse');

    title('DISTRIBUTION OF PRESSURE COEFFICIENT ALONG THE CHORD','FontName','Times New Roman', 'FontWeight','normal');
    cp_bot = min(cp);
    text(0.5, 0.90*cp_bot, sprintf('\\alpha = %g deg.', rad2deg(alpha)),'FontSize',11, 'FontName','Times New Roman');
    text(0.5, 0.75*cp_bot, sprintf('N_{div} = %d', Ndiv),'FontSize',11, 'FontName','Times New Roman');

end