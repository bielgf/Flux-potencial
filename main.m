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

computeGeometry2airfoil(16,4)