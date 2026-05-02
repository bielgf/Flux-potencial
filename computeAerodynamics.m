function [CL,L,CM1_4,M1_4,CM0,cp] = computeAerodynamics(Ndiv,geom,gamma,Qinfmod,rho)

    l     = geom.l;
    Xc    = geom.Xc;
    c     = geom.c;
    delta = geom.delta;
    Nc    = geom.Nc;

    cl    = zeros(Ndiv,1);
    cp    = zeros(Ndiv,1);
    cm1_4 = zeros(Ndiv,1);
    cm0   = zeros(Ndiv,1);

    for ii=1:Ndiv
        cl(ii)       = 2*gamma(ii)*l(ii)/Qinfmod;
        cp(ii)       = 1 - (gamma(ii)/Qinfmod)^2;
        cm1_4(ii)    = cp(ii)*((Xc(ii,1)/c - 0.25)*(delta(ii,1)/c) + (Xc(ii,2)/c)*(delta(ii,2)/c));
        cm0(ii)      = cp(ii)*((Xc(ii,1)/c)*(delta(ii,1)/c) + (Xc(ii,2)/c)*(delta(ii,2)/c));

        cp_comp(ii)  = cp(ii)/(sqrt(1 - Minf^2) + (Minf^2/(1 + sqrt(1 - Minf^2)))*(cp(ii)/2));
        cl_compr(ii) = cp_comp(ii)*l(ii)*Nc(ii);
    end
    
    CL    = sum(cl);
    L     = CL*0.5*Qinfmod^2*2*c*rho;
    CM1_4 = sum(cm1_4);
    M1_4  = CM1_4*0.5*rho*Qinfmod^2*c^2;
    CM0   = sum(cm0);

end