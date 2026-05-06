function geom = computeGeometry(Ndiv)

    % Open and read the file
    filename = sprintf('HQ_300\\HQ300_%.0f.txt',Ndiv);
    data     = load(filename);
    
    % Extract variables
    geom.X = [data(:, 2) data(:, 3)];
    geom.c = 1;
    
    % Preallocate vectors
    geom.l     = zeros(Ndiv,1);
    geom.Xc    = zeros(Ndiv,2);
    geom.delta = zeros(Ndiv,2);
    geom.ca    = zeros(Ndiv,1);
    geom.sa    = zeros(Ndiv,1);
    geom.Nc    = zeros(Ndiv,2);
    geom.Tc    = zeros(Ndiv,2);
    
    for jj = 1:Ndiv
        geom.l(jj)       = sqrt((geom.X(jj,1) - geom.X(jj+1,1))^2 + (geom.X(jj,2) - geom.X(jj+1,2))^2);     % Panel's lenght
        geom.Xc(jj,:)    = (geom.X(jj,:) + geom.X(jj+1,:))/2;                                               % Panel's geometric center
        geom.delta(jj,:) = geom.X(jj+1,:) - geom.X(jj,:);                                                   % Increments in X and Z
        geom.ca(jj)      = (geom.X(jj+1,1) - geom.X(jj,1))/geom.l(jj);                                      % Cosinus function
        geom.sa(jj)      = (geom.X(jj,2) - geom.X(jj+1,2))/geom.l(jj);                                      % Sinus function
        geom.Nc(jj,:)    = [geom.sa(jj,1),geom.ca(jj,1)];                                                   % Normal vectors coordinates
        geom.Tc(jj,:)    = [geom.ca(jj,1),-geom.sa(jj,1)];                                                  % Tangent vector coordinates
    end
    
end