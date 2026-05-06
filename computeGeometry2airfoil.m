function geom = computeGeometry2airfoil(Ndiv, delta_e_deg)

c1 = 0.64; % chord main airfoil
c2 = 0.34; % chord second airfoil
d  = 0.02; % gap

% Read NACA 0012 file
filename = sprintf('NACA0012/NACA_10_N_%d.txt', Ndiv);
data = load(filename);
Xbase = [data(:,2) data(:,3)];

% AIRFOIL 1 (horizontal stabilizer) 

X1 = [c1*Xbase(:,1), c1*Xbase(:,2)];

% AIRFOIL 2 (elevator)

X2 = [c2*Xbase(:,1), c2*Xbase(:,2)]; % Scale

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

figure; hold on; axis equal; grid on
plot(X1(:,1),X1(:,2),'b.-')
plot(X2(:,1),X2(:,2),'r.-')
legend('Main airfoil','Elevator')
title(sprintf('Two-element geometry, \\delta_e = %.1f°',delta_e_deg))
xlabel('x'); ylabel('z');