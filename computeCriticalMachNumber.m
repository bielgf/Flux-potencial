function Mcr = computeCriticalMachNumber(cp,gam)

    % Kármán-Tsien Compressibility Correction & Critical Mach Number
    Cp0     = cp;
    Cp0_min = min(Cp0);
    
    % Kármán-Tsien correction
    KT = @(Cp0_val, M) Cp0_val./(sqrt(1 - M.^2) + (M.^2./(1 + sqrt(1 - M.^2))).*(Cp0_val/2));
    
    % Critical Cp (local Mach = 1, isentropic)
    Cp_crit = @(M) (2./(gam.*M.^2)).*(((2/(gam+1)).*(1 + (gam-1)/2.*M.^2)).^(gam/(gam-1)) - 1 );
    
    residual = @(M) KT(Cp0_min, M) - Cp_crit(M);
    
    % Bracket search
    M_scan   = linspace(0.02, 0.99, 5000);
    denom_KT = sqrt(1 - M_scan.^2) + (M_scan.^2./(1 + sqrt(1 - M_scan.^2))).*(Cp0_min/2);
    
    valid    = denom_KT > 0; % physical region only
    M_valid  = M_scan(valid);
    res_scan = arrayfun(residual,M_valid);
    
    % Find first sign change
    sc = find(diff(sign(res_scan)) ~= 0, 1, 'first');
    
    if isempty(sc)
        error('No critical Mach found. Cp0_min = %.4f — body may have no suction.', Cp0_min);
    end
    
    Mcr = fzero(residual, [M_valid(sc), M_valid(sc+1)]);

end