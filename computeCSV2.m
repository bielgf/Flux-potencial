function gamma = computeCSV2(Ndiv_main,Ndiv_second,geom,Qinf)

    Tc   = geom.Tc;
    Xc   = geom.Xc;
    X    = geom.X;
    l    = geom.l;
    ca   = geom.ca;
    sa   = geom.sa;

    [A11, A12, B1] = computeA(Ndiv_main, Tc(1:Ndiv_main,:), Xc(1:Ndiv_main,:), X(1:Ndiv_main+1,:), ...
        l(1:Ndiv_main), ca(1:Ndiv_main), sa(1:Ndiv_main), Qinf, X(Ndiv_main+2:end,:), l(Ndiv_main+1:end), ca(Ndiv_main+1:end), sa(Ndiv_main+1:end));
    [A22, A21, B2] = computeA(Ndiv_second, Tc(Ndiv_main+1:end,:), Xc(Ndiv_main+1:end,:), X(Ndiv_main+2:end,:), ...
        l(Ndiv_main+1:end), ca(Ndiv_main+1:end), sa(Ndiv_main+1:end), Qinf, X(1:Ndiv_main+1,:), l(1:Ndiv_main), ca(1:Ndiv_main), sa(1:Ndiv_main));    

    a = [A11, A12; A21, A22];
    % a = [A11, zeros(Ndiv,1); zeros(1,Ndiv), 1];
    % a = [1, zeros(1,Ndiv); zeros(Ndiv,1), A22];
    b = [B1; B2];
    % b = [B1; 0];
    % b = [0;B2];

    % --- Kutta 1 ---
    % gamma_1 + gamma_Ndiv = 0 
    idx_K1 = fix(Ndiv_main/4)+1; 
    a(idx_K1, :) = 0;           % Clear row
    a(idx_K1, 1) = 1;           % Coeficient primer panell perfil 1
    a(idx_K1, Ndiv_main) = 1;   % Coeficient últim panell perfil 1
    b(idx_K1) = 0;

    % --- Kutta 2 ---
    % gamma_Ndiv+1 + gamma_Ntotal = 0 
    idx_K2 = Ndiv_main + fix(Ndiv_second/4);
    a(idx_K2, :) = 0;                        
    a(idx_K2, Ndiv_main+1) = 1;              
    a(idx_K2, Ndiv_main + Ndiv_second) = 1;  
    b(idx_K2) = 0;

    gamma         = a\b;
    gamma(idx_K1) = 0.5*(gamma(idx_K1-1)+gamma(idx_K1+1));
    gamma(idx_K2) = 0.5*(gamma(idx_K2-1)+gamma(idx_K2+1));

end

