alpha_vec = deg2rad(0:2:8);

Cl_vec = zeros(length(alpha_vec),1);
Cm_vec = zeros(length(alpha_vec),1);

for k = 1:length(alpha_vec)
    [Cl_vec(k),Cm_vec(k)] = vortexSolver(alpha_vec(k));
end