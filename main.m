%------------------------------------------------------------------%
%%%%% Flux Potencial. Aerodinàmica Numèrica en Perfils i Ales %%%%%%
%%%%%%%%%%%%%%% Judith Bailen, Eulàlia Caballol %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Biel González, Júlia Soliva  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% AMVO - Course 2025-2026 %%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------%

% ---------- PART 1: THE CONSTANT STRENGTH VORTEX METHOD --------- %
% ---------- APPLIED TO SIMPLE AND TWO-ELEMENT AIRFOILS ---------- %

clear;clc;close all

%% HQ300 airfoil study

alpha = deg2rad(0:2:8);
Ndiv  = [16, 32, 64, 128, 256, 512];

nA = length(alpha);
nN = length(Ndiv);

CL_table   = zeros(nA,nN);
CM14_table = zeros(nA,nN);
L_table    = zeros(nA,nN);
M14_table  = zeros(nA,nN);

CM0      = zeros(1,nA);  
AC_table = zeros(1,nN);

rho     = 1.225;
gam     = 1.4;
Qinfmod = 1;
Minf    = 0;


for jj = 1:nN
    for ii = 1:nA
        Qinf  = Qinfmod*[cos(alpha(ii)),sin(alpha(ii))];
        geom  = computeGeometry(Ndiv(jj));
        gamma = computeCSV(Ndiv(jj),geom,Qinf);
        [CL,L,CM1_4,M1_4,CM0(ii),cp,~] = computeAerodynamics(Ndiv(jj),geom,gamma,Qinfmod,rho,alpha(ii),Minf);
        CL_table(ii,jj)   = CL;
        CM14_table(ii,jj) = CM1_4;
        L_table(ii,jj)    = L;
        M14_table(ii,jj)  = M1_4;
        
        if Ndiv(jj) == 32  && alpha(ii) == deg2rad(8); plotGammaDistribution(Ndiv(jj),geom.X,geom.Nc,gamma); end
        if Ndiv(jj) == 256 && alpha(ii) == deg2rad(8); plotCpChordDistribution(geom.X,cp,alpha(ii),Ndiv(jj));
                                                       plotCpDistribution(geom.X,geom.Nc,geom.Xc,cp,alpha(ii),Ndiv(jj),CL); end
    end
    dCM0_da = polyfit(alpha(:),CM0(:),1);
    dCL_da  = polyfit(alpha(:),CL_table(:,jj),1);
    AC_table(jj) = -dCM0_da(1)/dCL_da(1);
end

plotCLandCM14vsalpha(alpha,CL_table,CM14_table);
plotConvergence(Ndiv,CL_table,alpha);


Mcr = zeros(1,nA-2);

for ii = 1:nA-2
    Qinf  = Qinfmod*[cos(alpha(ii)),sin(alpha(ii))];
    geom  = computeGeometry(512);
    gamma = computeCSV(512,geom,Qinf);
    [~,~,~,~,~,cp,~] = computeAerodynamics(512,geom,gamma,Qinfmod,rho,alpha(ii),Minf);
    Mcr(ii) = computeCriticalMachNumber(cp,gam);
end


Mcr_val      = [0, Mcr(end) - 0.15, Mcr(end) - 0.10, Mcr(end) - 0.05, Mcr(end)]; % The first being 0 to compare with the results obtained previously
CL_Mcr_table = zeros(1,length(Mcr_val));

for ii = 1:length(Mcr_val)
    Qinf  = Qinfmod*[cos(alpha(3)),sin(alpha(3))];
    geom  = computeGeometry(512);
    gamma = computeCSV(512,geom,Qinf);
    [~,~,~,~,~,~,CL_comp] = computeAerodynamics(512,geom,gamma,Qinfmod,rho,alpha(3),Mcr_val(ii));
    CL_Mcr_table(ii) = CL_comp;
end

%% Two NACA 0012 airfoils tandem
clear;clc;

% Disclaimer: It is important for this code to mantain the same number of
% Ndiv_main and Ndiv_second, as well as the same number of alpha and
% delta_e combinations, for how the loops are defined.

Ndiv_main   = [16, 32, 64, 128, 256, 512];
Ndiv_second = [16, 32, 64, 128, 256, 512];
Ntotal      = Ndiv_main + Ndiv_second;
subsection  = 5;  % 4 or 5

switch subsection
    case 4
        alpha   = deg2rad(0:2:8);
        delta_e = deg2rad([0;0;0;0;0]);
    case 5
        alpha   = deg2rad([4;4;4;4;4]);
        delta_e = deg2rad(0:4:16);
    otherwise
        disp('Incorrect exercise subsection')
end

nA = length(alpha);
nN = length(Ndiv_main);

CL_table   = zeros(nA,nN);
CM14_table = zeros(nA,nN);
L_table    = zeros(nA,nN);
M14_table  = zeros(nA,nN);

rho     = 1.225;
gam     = 1.4;
Qinfmod = 1;
Minf    = 0;

for jj = 1:nN
    for ii = 1:nA
        Qinf  = Qinfmod*[cos(alpha(ii)),sin(alpha(ii))];
        geom  = computeGeometry2airfoil(Ndiv_main(jj), Ndiv_second(jj), delta_e(ii));

        gamma = computeCSV2(Ndiv_main(jj),Ndiv_second(jj),geom,Qinf);
        [CL,L,CM1_4,M1_4,CM0,cp,~] = computeAerodynamics2(Ndiv_main(jj),Ndiv_second(jj),geom,gamma,Qinfmod,rho,alpha(ii),Minf);
        CL_table(ii,jj)   = CL;
        CM14_table(ii,jj) = CM1_4;
        L_table(ii,jj)    = L;
        M14_table(ii,jj)  = M1_4;     
    end

    if Ndiv_main(jj) == 256 && subsection == 4; plotCpChordDistribution2(geom.X,cp,alpha(ii),delta_e(ii),Ndiv_main(jj),Ndiv_second(jj));
                                                plotCpDistribution2(geom.X,geom.Nc,geom.Xc,cp,alpha(ii),delta_e(ii),Ndiv_main(jj),Ndiv_second(jj)); end
    if Ndiv_main(jj) == 256 && subsection == 5; plotCpChordDistribution2(geom.X,cp,alpha(ii),delta_e(ii),Ndiv_main(jj),Ndiv_second(jj));
                                                plotCpDistribution2(geom.X,geom.Nc,geom.Xc,cp,alpha(ii),delta_e(ii),Ndiv_main(jj),Ndiv_second(jj)); end 
end

% PLOTS
if subsection == 4; plotCLandCM14vsalpha(alpha,CL_table,CM14_table); end
if subsection == 5; plotCLandCM14vsdelta(delta_e,CL_table,CM14_table); end


%% 
% ---------- PART 2: PRANDTL’S LIFTING LINE MODEL --------------------- %
% ---------- APPLIED TO COMPOUND WINGS OF LARGE ASPECT RATIO ---------- %

%% Study of the wing isolated - HQ300

clear;clc;close all

% Input Data

% Geometric Data
b    = 15;
b_h  = 3;
c_r  = 0.95;
c_rh = 0.5;
c_t  = 0.55;
c_th = 0.3;
l_h  = 4;
l_v  = 1.2;
Sv   = 1.5;
i_w  = 0;
i_h  = deg2rad(3);
Sw   = (c_r + c_t)/2*b;

% Numerical Data
Nw = 512;
Nh = Nw/4;

% Aerodynamic Data
rho       = 1.225;
alpha     = deg2rad(4);
theta_tip = deg2rad(0);
Claw      = 6.935845305522705;  
Cl0w      = 0.623151006958811;
i_inf     = [cos(alpha) , 0 , sin(alpha)];
ur        = -i_inf;
k_inf     = [-sin(alpha) , 0 , cos(alpha)];
Qinf_mod  = norm(i_inf);

% Geometric discretization of the wing-canard configuration
dy_w = b/Nw;
dy_h = b_h/Nh;
yP_w = -b/2:dy_w:b/2;
yP_h = -b_h/2:dy_h:b_h/2;

P_w     = [zeros(size(yP_w)) ; yP_w ; zeros(size(yP_w))]';
P_h     = [l_h*ones(size(yP_h)) ; yP_h ; zeros(size(yP_h))]';
P_w_mid = [zeros(1,size(yP_w,2)-1) ; (yP_w(1:end-1)+yP_w(2:end))/2 ; zeros(1,size(yP_w,2)-1)]';
P_h_mid = [l_h*ones(1,size(yP_h,2)-1) ; (yP_h(1:end-1)+yP_h(2:end))/2 ; zeros(1,size(yP_h,2)-1)]';

theta_mid = theta_tip*(2*abs(P_w_mid(:,2))/b);
cwi05     = c_r + (c_t - c_r)*(2*abs(P_w_mid(:,2))/b);
chi05     = c_rh + (c_th - c_rh)*(2*abs(P_h_mid(:,2))/b_h);



gamma = computeWingConfig(Nw,P_w,P_w_mid,ur,cwi05,k_inf,Claw,Qinf_mod,Cl0w,alpha,theta_mid);

L  = rho*Qinf_mod*sum(gamma*dy_w);
CL = L/(0.5*rho*Qinf_mod^2*Sw);


function gamma = computeWingConfig(Nw,P_w,P_w_mid,ur,cwi05,k_inf,Cla,Q_inf_mod,Cl0,alpha,thetai05) % Punt 1 de la part 2
    
    a = zeros(Nw,Nw);
    b = zeros(Nw,1);

    for ii = 1:Nw
        b(ii,1) = 0.5*cwi05(ii)*Q_inf_mod*(Cl0 + Cla*(alpha + thetai05(ii)));
        for jj = 1:Nw
            if ii == jj
                v = computeHSV(ii,jj,P_w(jj,:),P_w(jj+1,:),P_w_mid(ii,:),ur,'WingWing');
                a(ii,jj) = -0.5*Cla*cwi05(ii)*dot(v,k_inf) + 1;
            else
                v = computeHSV(ii,jj,P_w(jj,:),P_w(jj+1,:),P_w_mid(ii,:),ur,'WingWing');
                a(ii,jj) = -0.5*Cla*cwi05(ii)*dot(v,k_inf);
            end
        end
    end

    gamma = a\b;

end


function gamma = computeWingCanardConfig(Nw,Nh)

    a = zeros(Nw+Nh,Nw+Nh);
    b = zeros(Nw+Nh,1);

    for ii = 1:Nw
        b(ii,1) = 0.5*cwi05(ii)*Q_inf_mod*(Cl0w + Claw*(alpha + thetai05(ii)));
        for jj = 1:Nw
            if ii == jj
                v = computeHSV(ii,jj,P_w(jj,:),P_w(jj+1,:),P_w_mid(ii,:),ur,'WingWing');
                a(ii,jj) = -0.5*Claw*cwi05(ii)*dot(v,k_inf) + 1;
            else
                v = computeHSV(ii,jj,P_w(jj,:),P_w(jj+1,:),P_w_mid(ii,:),ur,'WingWing');
                a(ii,jj) = -0.5*Claw*cwi05(ii)*dot(v,k_inf);
            end
        end

        for jj = 1:Nh
            v = computeHSV(ii,jj,P_h(jj,:),P_h(jj+1,:),P_w_mid(ii,:),ur,'WingCnrd');
            a(ii,Nw+jj) = -0.5*Claw*cwi05(ii)*dot(v,k_inf);
        end
    end

    for ii = 1:Nh
        b(Nw+ii,1) = 0.5*chi05(ii)*Q_inf_mod*(Cl0h + Clah*alpha + i_h);
        for jj = 1:Nh
            if ii == jj
                v = computeHSV(ii,jj,P_h(jj,:),P_h(jj+1,:),P_h_mid(ii,:),ur,'CnrdCnrd');
                a(Nw+ii,Nw+jj) = -0.5*Clah*chi05(ii)*dot(v,k_inf) + 1;
            else
                v = computeHSV(ii,jj,P_h(jj,:),P_h(jj+1,:),P_h_mid(ii,:),ur,'CnrdCnrd');
                a(Nw+ii,Nw+jj) = -0.5*Clah*chi05(ii)*dot(v,k_inf);
            end
        end

        for jj = 1:Nw
            v = computeHSV(ii,jj,P_w(jj,:),P_w(jj+1,:),P_h_mid(ii,:),ur,'CnrdWing');
            a(Nw+ii,jj) = -0.5*Clah*chi05(ii)*dot(v,k_inf);
        end
    end

    gamma = a\b;

end