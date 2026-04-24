%% Preprocess

Ndiv   = 500;
deltaB = 2*pi/Ndiv;
R      = 0.5;

Nc = zeros(Ndiv,2);
Tc = zeros(Ndiv,2);

ca = zeros(Ndiv,1);
sa = zeros(Ndiv,1);
l  = zeros(Ndiv,1);

X  = zeros(Ndiv+1,2);
Xc = zeros(Ndiv,2);

cl = zeros(Ndiv,1);
cp = zeros(Ndiv,1);

for ii = 1:Ndiv+1
    X(ii,:)  = R*[cos((ii-1)*deltaB), -sin((ii-1)*deltaB)];
end

for jj = 1:Ndiv
    l(jj)    = sqrt((X(jj,1) - X(jj+1,1))^2 + (X(jj,2) - X(jj+1,2))^2);
    Xc(jj,:) = (X(jj,:) + X(jj+1,:))/2;
    ca(jj)   = (X(jj+1,1) - X(jj,1))/l(jj);
    sa(jj)   = (X(jj,2) - X(jj+1,2))/l(jj);

    Nc(jj,:) = [sa(jj,1),ca(jj,1)];
    Tc(jj,:) = [ca(jj,1),-sa(jj,1)];
end

%% Process

a = zeros(Ndiv,Ndiv);
b = zeros(Ndiv,1);

alpha = deg2rad(0);

rho = 1.225;

Qinfmod = 1;
Qinf    = Qinfmod*[cos(alpha),sin(alpha)];

for ii = 1:Ndiv
    b(ii) = -dot(Qinf,Tc(ii,:));
    for jj = 1:Ndiv
        if jj ~= ii
            Xc_ij_pan = (Xc(ii,1) - X(jj,1))*ca(jj) - (Xc(ii,2) - X(jj,2))*sa(jj);
            Zc_ij_pan = (Xc(ii,1) - X(jj,1))*sa(jj) + (Xc(ii,2) - X(jj,2))*ca(jj);

            r1 = sqrt(Xc_ij_pan^2 + Zc_ij_pan^2);
            r2 = sqrt((Xc_ij_pan - l(jj))^2 + Zc_ij_pan^2);

            theta1 = atan2(Zc_ij_pan,Xc_ij_pan);
            theta2 = atan2(Zc_ij_pan,Xc_ij_pan - l(jj));

            u_ij_pan = (theta2 - theta1)/(2*pi);
            w_ij_pan = (1/(4*pi))*log(r2^2/r1^2);

            u_ij = u_ij_pan*ca(jj) + w_ij_pan*sa(jj);
            w_ij = -u_ij_pan*sa(jj) + w_ij_pan*ca(jj);

            a(ii,jj) = dot([u_ij,w_ij],Tc(ii,:));
        else
            a(ii,jj) = -0.5;
        end
    end
end

K         = fix(Ndiv/4);
a(K,:)    = 0;
a(K,1)    = 1;
a(K,Ndiv) = 1;
b(K)      = 0;

gamma = a\b;
gamma(K) = 0.5*(gamma(K-1)+gamma(K+1));

for ii=1:Ndiv
     cl(ii) = 2*gamma(ii)*l(ii)/Qinfmod;
     cp(ii) = 1 - (gamma(ii)/Qinfmod)^2;
end

CL = sum(cl);
L  = CL*0.5*Qinfmod^2*2*R*rho;

theta_c = atan2(Xc(:,2), Xc(:,1));
[theta_sorted, idx] = sort(theta_c);
figure;
plot(rad2deg(theta_sorted), cp(idx), '-o','MarkerSize',3);
xlabel('\theta (rad)'); ylabel('C_p');
title('Distribución de C_p en los paneles (ordenada por \theta)');
grid on;


% Pel segon exercici part del treball - Part 1
% A = [A11 A12 ; A21 A22]

% Hi han dues condicions de Kutta, ens hem de carregar dues files de la
% matriu A (corresponent a panells del intradós dels perfils corresponents)