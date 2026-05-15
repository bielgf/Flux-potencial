function [CL,L,CM1_4,M1_4,CM0,cp,CL_comp] = computeAerodynamics2(Ndiv1,Ndiv2,geom,gamma,Qinfmod,rho,alpha,Minf)

    Ntotal = Ndiv1 + Ndiv2;

    l     = geom.l;
    Xc    = geom.Xc;
    c     = geom.c;
    delta = geom.delta;
    Nc    = geom.Nc;

    cl      = zeros(Ntotal,1);
    cp      = zeros(Ntotal,1);
    cm1_4   = zeros(Ntotal,1);
    cm0     = zeros(Ntotal,1);
    cp_comp = zeros(Ntotal,1);
    cl_comp = zeros(Ntotal,2);

    for ii=1:Ntotal
        cl(ii)       = 2*gamma(ii)*l(ii)/Qinfmod;
        cp(ii)       = 1 - (gamma(ii)/Qinfmod)^2;
        cm1_4(ii)    = cp(ii)*((Xc(ii,1)/c - 0.25)*(delta(ii,1)/c) + (Xc(ii,2)/c)*(delta(ii,2)/c));
        cm0(ii)      = cp(ii)*((Xc(ii,1)/c)*(delta(ii,1)/c) + (Xc(ii,2)/c)*(delta(ii,2)/c));

        cp_comp(ii)    = cp(ii)/(sqrt(1 - Minf^2) + (Minf^2/(1 + sqrt(1 - Minf^2)))*(cp(ii)/2));
        cl_comp(ii,1) = cp_comp(ii)*l(ii)*Nc(ii,1);
        cl_comp(ii,2) = cp_comp(ii)*l(ii)*Nc(ii,2);
    end
    
    CL    = sum(cl);
    L     = CL*0.5*Qinfmod^2*2*c*rho;
    CM1_4 = sum(cm1_4);
    M1_4  = CM1_4*0.5*rho*Qinfmod^2*c^2;
    CM0   = sum(cm0);

    CL_comp = sum(cl_comp);
    k_inf   = [-sin(alpha), cos(alpha)];                % Perpendicular vector to freestream velocity
    CL_comp = -(k_inf(1,1)*CL_comp(1,1) + k_inf(1,2)*CL_comp(1,2));

end
