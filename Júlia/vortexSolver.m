function [Cl_total,Cm1_4_total] = vortexSolver(alpha)

% PRE-PROCÉS

coords = load('HQ300_16.txt');   % fitxer amb punts del perfil
X = coords; % punts nodals [x z]

% Assegurar que el perfil està tancat
if X(1,:) ~= X(end,:)
    X(end+1,:) = X(1,:);
end

N = size(X,1)-1; % nombre de panells
Xc = zeros(N,2); % punts de control
l  = zeros(N,1); % longituds dels panells
delta = zeros(N,2); % increments

ca = zeros(N,1);
sa = zeros(N,1);
nc = zeros(N,2);
tc = zeros(N,2);

for j = 1:N
    Xc(j,:) = (X(j,:) + X(j+1,:))/2; % coordenades x i z dels punts de control

    l(j) = sqrt((X(j,1) - X(j+1,1))^2 + (X(j,2) - X(j+1,2))^2);
    delta(j,:) = X(j+1,:) - X(j,:);

    ca(j) = (X(j+1,1)-X(j,1))/l(j); % cosinus de alpha
    sa(j) = (X(j,2)-X(j+1,2))/l(j); % sinus de alpha
    nc(j,:) = [sa(j), ca(j)]; % vector normal
    tc(j,:) = [ca(j), -sa(j)]; % vector tangencial
end

% PROCÉS

a = zeros(N,N);
b = zeros(N,1);
Qinf = 1*[cos(alpha), sin(alpha)]; % el valor de 1 és arbitrari

for j = 1:N
    b(j) = -dot(Qinf,tc(j,:));
end

for i = 1:N
    for j = 1:N
        if j ~= i
            Xcij_pan = (Xc(i,1) - X(j,1))*ca(j) - (Xc(i,2) - X(j,2))*sa(j); 
            Zcij_pan = (Xc(i,1) - X(j,1))*sa(j) + (Xc(i,2) - X(j,2))*ca(j);
            r1 = sqrt(Xcij_pan^2 + Zcij_pan^2);
            r2 = sqrt((Xcij_pan-l(j))^2 + Zcij_pan^2);
            theta1 = atan2(Zcij_pan,Xcij_pan);
            theta2 = atan2(Zcij_pan, Xcij_pan-l(j));
            Uij_pan = log((r2^2)/(r1^2))/(4*pi);
            Wij_pan = (theta2 - theta1) / (2*pi);
            U(i,j) = (Uij_pan*ca(j) + Wij_pan*sa(j));
            W(i,j) = (-Uij_pan*sa(j) + Wij_pan*ca(j));
            a(i,j) = dot([U(i,j),W(i,j)],tc(i,:));
        else
            a(i,j)= -0.5;
        end
    end
end

% gamma_nokutta = a\b;

% Kutta condition
a(N,:) = 0;
a(N,1) = 1;
a(N,N) = 1;
b(N)   = 0;

gamma = a\b; % vector columna

% POST-PROCÉS

rho = 1.225;
c = max(X(:,1)) - min(X(:,1)); % corda
gammatotal = gamma'*l; % sumatori de l(j) * gamma(j)

Cl  = zeros(N,1);
Cp  = zeros(N,1);
Cm1_4 = zeros(N,1);

for i = 1:N
    Cl(i)    = (2*gamma(i)*l(i))/(norm(Qinf)*c);
    Cp(i,1) = 1-(gamma(i)/norm(Qinf))^2;
    Cm1_4(i) = Cp(i)*((Xc(i,1)/c)*(delta(i,1)/c) + (Xc(i,2)/c)*(delta(i,2)/c)) - 0.25*Cl(i);
end

Cl_total = sum(Cl);
% Lift = 0.5*rho*norm(Qinf)^2*Cl_total*c;
Cm1_4_total = sum(Cm1_4);
% M1_4  = 0.5*Cm1_4*rho*norm(Qinf)^2*c;

end
