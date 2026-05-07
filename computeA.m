function [A11,A12,B1] = computeA(Ndiv, Tc1, Xc1, X1, l1, ca1, sa1, Qinf, X2, l2, ca2, sa2)

    A11  = zeros(Ndiv,Ndiv);
    A12  = zeros(Ndiv,Ndiv);
    B1   = zeros(Ndiv,1);

    for ii = 1:Ndiv
        B1(ii) = -dot(Qinf,Tc1(ii,:));

        for jj = 1:Ndiv % Aii matrix
            if jj ~= ii
                Xc_ij_pan = (Xc1(ii,1) - X1(jj,1))*ca1(jj) - (Xc1(ii,2) - X1(jj,2))*sa1(jj);
                Zc_ij_pan = (Xc1(ii,1) - X1(jj,1))*sa1(jj) + (Xc1(ii,2) - X1(jj,2))*ca1(jj);

                r1 = sqrt(Xc_ij_pan^2 + Zc_ij_pan^2);
                r2 = sqrt((Xc_ij_pan - l1(jj))^2 + Zc_ij_pan^2);

                theta1 = atan2(Zc_ij_pan,Xc_ij_pan);
                theta2 = atan2(Zc_ij_pan,Xc_ij_pan - l1(jj));

                u_ij_pan = (theta2 - theta1)/(2*pi);
                w_ij_pan = (1/(4*pi))*log(r2^2/r1^2);

                u_ij = u_ij_pan*ca1(jj) + w_ij_pan*sa1(jj);
                w_ij = -u_ij_pan*sa1(jj) + w_ij_pan*ca1(jj);

                A11(ii,jj) = dot([u_ij,w_ij],Tc1(ii,:));
            else
                A11(ii,jj) = -0.5;
            end
        end

        for jj = 1:Ndiv % Aij matrix
            if jj ~= ii
                Xc_ij_pan = (Xc1(ii,1) - X2(jj,1))*ca2(jj) - (Xc1(ii,2) - X2(jj,2))*sa2(jj);
                Zc_ij_pan = (Xc1(ii,1) - X2(jj,1))*sa2(jj) + (Xc1(ii,2) - X2(jj,2))*ca2(jj);

                r1 = sqrt(Xc_ij_pan^2 + Zc_ij_pan^2);
                r2 = sqrt((Xc_ij_pan - l2(jj))^2 + Zc_ij_pan^2);

                theta1 = atan2(Zc_ij_pan,Xc_ij_pan);
                theta2 = atan2(Zc_ij_pan,Xc_ij_pan - l2(jj));

                u_ij_pan = (theta2 - theta1)/(2*pi);
                w_ij_pan = (1/(4*pi))*log(r2^2/r1^2);

                u_ij = u_ij_pan*ca2(jj) + w_ij_pan*sa2(jj);
                w_ij = -u_ij_pan*sa2(jj) + w_ij_pan*ca2(jj);

                A12(ii,jj) = dot([u_ij,w_ij],Tc1(ii,:));
            end
        end

    end

end

