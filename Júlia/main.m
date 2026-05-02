alpha_vec = deg2rad(0:2:8); % degrees or radiants

Cl_vec = zeros(length(alpha_vec),1);
Cm_vec = zeros(length(alpha_vec),1);
Mcr = zeros(length(alpha_vec),1);

for k = 1:length(alpha_vec)
    [Cl_vec(k),Cm_vec(k),Mcr(k)] = vortexSolver(alpha_vec(k));
end

alpha_4rad = 4*(pi/180);
Cl0_4deg = Cl_vec(alpha_vec==alpha_4rad);
Mcr_4deg = Mcr(alpha_vec==alpha_4rad);
M_list = [Mcr_4deg-0.15, Mcr_4deg-0.10, Mcr_4deg-0.05, Mcr_4deg ];

gamma_air = 1.4;
CL_comp = zeros(length(M_list),1);

for k = 1:length(M_list)

    M = M_list(k);

    CL_comp(k) = Cl0_4deg / (sqrt(1-M^2)+(Cl0_4deg/2)*(M^2/(1+sqrt(1-M^2))));

end

disp('Lift coefficient compressible per alpha = 4º')
disp(table(M_list', CL_comp,'VariableNames',{'Mach','Cl_compressible'}))