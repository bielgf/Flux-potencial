function plotCpDistribution2(X, Nc, Xc, cp, alpha, delta, Ndiv1, Ndiv2)

    delta = deg2rad(delta);

    % --- AJUSTE DE DIMENSIONES ---
    cp = cp(:);
    totalPanels = Ndiv1 + Ndiv2;
    if length(cp) > totalPanels, cp = cp(1:totalPanels); end

    % --- CÁLCULO DEL OFFSET (Alineación vertical) ---
    % Borde de salida perfil 1 (X máximo del primer bloque)
    x_te1 = max(X(1:Ndiv1+1, 1));
    % Borde de ataque perfil 2 (X mínimo del segundo bloque)
    x_le2 = min(X(Ndiv1+2:end, 1));
    % Desplazamiento necesario
    dx_offset = x_te1 - x_le2;

    % Aplicamos el desplazamiento a las coordenadas originales
    X_plot = X; Xc_plot = Xc; Nc_plot = Nc;
    idx_nodes2 = Ndiv1+2 : size(X,1);
    idx_panels2 = Ndiv1+1 : totalPanels;
    
    X_plot(idx_nodes2, 1) = X_plot(idx_nodes2, 1) + dx_offset;
    Xc_plot(idx_panels2, 1) = Xc_plot(idx_panels2, 1) + dx_offset;

    % --- ROTACIÓN DEL SEGUNDO PERFIL ---
    % Punto de charnela: Ahora es el LE2 ya desplazado
    x_h = min(X_plot(idx_nodes2, 1));
    y_h = X_plot(Ndiv1+2, 2); % O el punto que prefieras como eje
    
    if delta ~= 0
        R = [cos(-delta), -sin(-delta); sin(-delta), cos(-delta)];
        
        % Rotar Nodos
        X_sub = (X_plot(idx_nodes2, :) - [x_h, y_h]) * R';
        X_plot(idx_nodes2, :) = X_sub + [x_h, y_h];
        
        % Rotar Centros
        Xc_sub = (Xc_plot(idx_panels2, :) - [x_h, y_h]) * R';
        Xc_plot(idx_panels2, :) = Xc_sub + [x_h, y_h];
        
        % Rotar Normales
        Nc_plot(idx_panels2, :) = Nc_plot(idx_panels2, :) * R';
    end

    % --- DIBUJO ---
    figure('Name','Cp Aligned','Color','w'); hold on; axis equal off;
    plot(X_plot(1:Ndiv1+1, 1), X_plot(1:Ndiv1+1, 2), 'k-', 'LineWidth', 1.5);
    plot(X_plot(idx_nodes2, 1), X_plot(idx_nodes2, 2), 'k-', 'LineWidth', 1.5);
    
    sc4 = 0.2 / max(abs(cp)); 
    mask_suc = cp < 0; mask_pres = cp >= 0;
    
    if any(mask_suc)
        quiver(Xc_plot(mask_suc,1), Xc_plot(mask_suc,2), ...
               abs(cp(mask_suc)).*sc4.*Nc_plot(mask_suc,1), ...
               abs(cp(mask_suc)).*sc4.*Nc_plot(mask_suc,2), 0, 'r');
    end
    if any(mask_pres)
        x0_b = Xc_plot(mask_pres,1) + cp(mask_pres).*sc4.*Nc_plot(mask_pres,1);
        y0_b = Xc_plot(mask_pres,2) + cp(mask_pres).*sc4.*Nc_plot(mask_pres,2);
        quiver(x0_b, y0_b, Xc_plot(mask_pres,1)-x0_b, Xc_plot(mask_pres,2)-y0_b, 0, 'b');
    end
    title('PRESSURE DISTRIBUTION');

    text(min(X(:,1))+0.5, max(X(:,2))+0.2, sprintf('\\alpha=%gº, \\delta=%gº', rad2deg(alpha), rad2deg(delta)));
end
