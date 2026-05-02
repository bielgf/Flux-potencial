function gamma = computeCSV(Ndiv,geom,Qinf)

    a  = zeros(Ndiv,Ndiv);
    b  = zeros(Ndiv,1);
    Tc = geom.Tc;
    Xc = geom.Xc;
    X  = geom.X;
    l  = geom.l;
    ca = geom.ca;
    sa = geom.sa;

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
    
    gamma    = a\b;
    gamma(K) = 0.5*(gamma(K-1)+gamma(K+1));

end