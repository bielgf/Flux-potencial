function gamma = computeCSV2(Ndiv,geom,Qinf)

    Ntotal = 2*Ndiv;

    Tc   = geom.Tc;
    Xc   = geom.Xc;
    X    = geom.X;
    l    = geom.l;
    ca   = geom.ca;
    sa   = geom.sa;

    [A11, A12, B1] = computeA(Ndiv, Tc(1:Ndiv,:), Xc(1:Ndiv,:), X(1:Ndiv,:), ...
        l(1:Ndiv), ca(1:Ndiv), sa(1:Ndiv), Qinf, X(Ndiv+1:end,:), l(Ndiv+1:end), ca(Ndiv+1:end), sa(Ndiv+1:end));
    [A22, A21, B2] = computeA(Ndiv, Tc(Ndiv+1:end,:), Xc(Ndiv+1:end,:), X(Ndiv+1:end,:), ...
        l(Ndiv+1:end), ca(Ndiv+1:end), sa(Ndiv+1:end), Qinf, X(1:Ndiv,:), l(1:Ndiv), ca(1:Ndiv), sa(1:Ndiv));    

    a = [A11, A12; A21, A22];
    b = [B1; B2];

    % --- Kutta 1 ---
    % gamma_1 + gamma_Ndiv = 0 
    idx_K1 = fix(Ndiv/4); 
    a(idx_K1, :) = 0;      % Clear row
    a(idx_K1, 1) = 1;      % Coeficient primer panell perfil 1
    a(idx_K1, Ndiv) = 1;   % Coeficient últim panell perfil 1
    b(idx_K1) = 0;

    % --- Kutta 2 ---
    % gamma_Ndiv+1 + gamma_Ntotal = 0 
    idx_K2 = Ndiv + fix(Ndiv/4);
    a(idx_K2, :) = 0;      % Clear row
    a(idx_K2, Ndiv+1) = 1; % Coeficient primer panell perfil 1
    a(idx_K2, Ntotal) = 1; % Coeficient últim panell perfil 1
    b(idx_K2) = 0;

    gamma         = a\b;
    gamma(idx_K1) = 0.5*(gamma(idx_K1-1)+gamma(idx_K1+1));
    gamma(idx_K2) = 0.5*(gamma(idx_K2-1)+gamma(idx_K2+1));

end

