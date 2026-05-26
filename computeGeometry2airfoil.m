function geom = computeGeometry2airfoil(Ndiv_main, Ndiv_second, delta_e_deg)

c1 = 0.64; % chord main airfoil
c2 = 0.34; % chord second airfoil
d  = 0.02; % gap

% MAIN AIRFOIL: Read NACA 0012 file
filename_main = sprintf('NACA0012\\NACA_10_N_%d.txt', Ndiv_main);
data_main = load(filename_main);
Xbase_main = [data_main(:,2) data_main(:,3)];

% SECOND AIRFOIL: Read NACA 0012 file
filename_second = sprintf('NACA0012\\NACA_10_N_%d.txt', Ndiv_second);
data_second = load(filename_second);
Xbase_second = [data_second(:,2) data_second(:,3)];

% MAIN AIRFOIL (horizontal stabilizer) 

X1 = [c1*Xbase_main(:,1), c1*Xbase_main(:,2)];

% SECOND AIRFOIL (elevator)

X2 = [c2*Xbase_second(:,1), c2*Xbase_second(:,2)]; % Scale

X2(:,1) = X2(:,1) + c1 + d;      % Translate d distance

delta_e = deg2rad(delta_e_deg);  % Deflection in radians 

xLE = c1 + d;
zLE = 0;

for i = 1:size(X2,1)     % Rotation
    x = X2(i,1) - xLE;
    z = X2(i,2) - zLE;
    
    X2(i,1) =  x*cos(delta_e) + z*sin(delta_e) + xLE;
    X2(i,2) = -x*sin(delta_e) + z*cos(delta_e) + zLE;
end

% Airfoils altogether

geom.X = [X1; X2];
geom.c = c1 + c2 + d;

Ntotal = Ndiv_main + Ndiv_second; 
geom.l     = zeros(Ntotal,1);
geom.Xc    = zeros(Ntotal,2);
geom.delta = zeros(Ntotal,2);
geom.ca    = zeros(Ntotal,1);
geom.sa    = zeros(Ntotal,1);
geom.Nc    = zeros(Ntotal,2);
geom.Tc    = zeros(Ntotal,2);

for jj = 1:Ndiv_main

    % Main airfoil (horizontal stabilizer)
    geom.l(jj)       = sqrt((X1(jj,1) - X1(jj+1,1))^2 + (X1(jj,2) - X1(jj+1,2))^2);          % Panel's lenght
    geom.Xc(jj,:)    = (X1(jj,:) + X1(jj+1,:))/2;                                            % Panel's geometric center
    geom.delta(jj,:) = X1(jj+1,:) - X1(jj,:);                                                % Increments in X and Z
    geom.ca(jj)      = (X1(jj+1,1) - X1(jj,1))/geom.l(jj);                                   % Cosinus function
    geom.sa(jj)      = (X1(jj,2) - X1(jj+1,2))/geom.l(jj);                                   % Sinus function
    geom.Nc(jj,:)    = [geom.sa(jj,1),geom.ca(jj,1)];                                        % Normal vectors coordinates
    geom.Tc(jj,:)    = [geom.ca(jj,1),-geom.sa(jj,1)];                                       % Tangent vector coordinates

end

for jj = 1:Ndiv_second

    % Second airfoil (elevator)
    geom.l(jj+Ndiv_main)       = sqrt((X2(jj,1) - X2(jj+1,1))^2 + (X2(jj,2) - X2(jj+1,2))^2);     % Panel's lenght
    geom.Xc(jj+Ndiv_main,:)    = (X2(jj,:) + X2(jj+1,:))/2;                                       % Panel's geometric center
    geom.delta(jj+Ndiv_main,:) = X2(jj+1,:) - X2(jj,:);                                           % Increments in X and Z
    geom.ca(jj+Ndiv_main)      = (X2(jj+1,1) - X2(jj,1))/geom.l(jj+Ndiv_main);                         % Cosinus function
    geom.sa(jj+Ndiv_main)      = (X2(jj,2) - X2(jj+1,2))/geom.l(jj+Ndiv_main);                         % Sinus function
    geom.Nc(jj+Ndiv_main,:)    = [geom.sa(jj+Ndiv_main,1),geom.ca(jj+Ndiv_main,1)];                         % Normal vectors coordinates
    geom.Tc(jj+Ndiv_main,:)    = [geom.ca(jj+Ndiv_main,1),-geom.sa(jj+Ndiv_main,1)];                        % Tangent vector coordinates
    
end

figure; hold on; axis equal; grid on
plot(X1(:,1),X1(:,2),'b.-')
plot(X2(:,1),X2(:,2),'r.-')
legend('Main airfoil','Elevator')
title(sprintf('Two-element geometry, \\delta_e = %.1f°',delta_e_deg))
xlabel('x'); ylabel('z');
