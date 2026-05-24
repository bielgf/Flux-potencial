function plotCpChordDistribution2(X, cp, alpha, delta, Ndiv1, Ndiv2)

    % Forzamos que cp sea un vector columna para evitar errores de dimensiones
    cp = cp(:); 

    % --- Segmentación de Datos ---
    % Perfil 1: Nodos 1 a Ndiv1+1. CP 1 a Ndiv1
    X1 = X(1:Ndiv1+1, 1);
    cp1 = cp(1:Ndiv1);
    
    % Perfil 2: Nodos Ndiv1+2 a Final. CP empieza en Ndiv1+1
    % Usamos Ndiv1+Ndiv2 para asegurar que el tamaño sea EXACTAMENTE Ndiv2
    X2 = X(Ndiv1+2:Ndiv1+Ndiv2+2, 1);
    cp2 = cp(Ndiv1+1:Ndiv1+Ndiv2);

    % Centros de paneles (X_mid1 tendrá longitud Ndiv1, X_mid2 longitud Ndiv2)
    X_mid1 = (X1(1:end-1) + X1(2:end)) / 2;
    X_mid2 = (X2(1:end-1) + X2(2:end)) / 2;

    % --- Lógica de Intradós/Extradós ---
    [~, LE_idx1] = min(X_mid1);
    idx_int1 = 1:LE_idx1;
    idx_ext1 = LE_idx1+1:length(X_mid1);

    [~, LE_idx2] = min(X_mid2);
    idx_int2 = 1:LE_idx2;
    idx_ext2 = LE_idx2+1:length(X_mid2);

    % --- Gráfica ---
    figure; hold on; grid on;
    
    % Perfil 1
    plot(X_mid1, cp1, 'k-', 'HandleVisibility', 'off');
    p1_int = plot(X_mid1(idx_int1), cp1(idx_int1), 'bo', 'MarkerSize', 4);
    p1_ext = plot(X_mid1(idx_ext1), cp1(idx_ext1), 'ro', 'MarkerSize', 4);
    
    % Perfil 2 (Línea discontinua)
    plot(X_mid2, cp2, 'k--', 'HandleVisibility', 'off');
    p2_int = plot(X_mid2(idx_int2), cp2(idx_int2), 'bs', 'MarkerSize', 4);
    p2_ext = plot(X_mid2(idx_ext2), cp2(idx_ext2), 'rs', 'MarkerSize', 4);

    xlabel('X', 'FontSize', 12);
    ylabel('C_p', 'FontSize', 12);
    legend([p1_int, p1_ext, p2_int, p2_ext], ...
           {'Int. P1', 'Ext. P1', 'Int. P2', 'Ext. P2'}, 'Location', 'northeast');
    
    set(gca, 'YDir', 'reverse');
    title('PRESSURE COEFFICIENT DISTRIBUTION','FontName','Times New Roman');
    
    cp_min = min(cp);
    text(0.1, 0.90*cp_min, sprintf('\\alpha = %gº | \\delta = %gº ', rad2deg(alpha), rad2deg(delta)), 'FontName','Times New Roman');
end