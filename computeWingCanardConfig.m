function [gamma,gamma_w,gamma_h] = computeWingCanardConfig(Nw,Nh,P_w,P_h,P_w_mid,P_h_mid,cwi05,chi05,Q_inf_mod,Cl0w,Cl0h,Claw,Clah,thetai05,i_h,alpha,ur,k_inf)

    a = zeros(Nw+Nh,Nw+Nh);
    b = zeros(Nw+Nh,1);

    for ii = 1:Nw
        b(ii,1) = 0.5*cwi05(ii)*Q_inf_mod*(Cl0w + Claw*(alpha + thetai05(ii)));
        % Wing-Wing influence
        for jj = 1:Nw
            if ii == jj
                v = computeHSV(ii,jj,P_w(jj,:),P_w(jj+1,:),P_w_mid(ii,:),ur,'WingWing');
                a(ii,jj) = -0.5*Claw*cwi05(ii)*dot(v,k_inf) + 1;
            else
                v = computeHSV(ii,jj,P_w(jj,:),P_w(jj+1,:),P_w_mid(ii,:),ur,'WingWing');
                a(ii,jj) = -0.5*Claw*cwi05(ii)*dot(v,k_inf);
            end
        end

        % Wing-Canard influence
        for jj = 1:Nh
            v = computeHSV(ii,jj,P_h(jj,:),P_h(jj+1,:),P_w_mid(ii,:),ur,'WingCnrd');
            a(ii,Nw+jj) = -0.5*Claw*cwi05(ii)*dot(v,k_inf);
        end
    end

    for ii = 1:Nh
        b(Nw+ii,1) = 0.5*chi05(ii)*Q_inf_mod*(Cl0h + Clah*(alpha + i_h));
        % Canard-Canard influence
        for jj = 1:Nh
            if ii == jj
                v = computeHSV(ii,jj,P_h(jj,:),P_h(jj+1,:),P_h_mid(ii,:),ur,'CnrdCnrd');
                a(Nw+ii,Nw+jj) = -0.5*Clah*chi05(ii)*dot(v,k_inf) + 1;
            else
                v = computeHSV(ii,jj,P_h(jj,:),P_h(jj+1,:),P_h_mid(ii,:),ur,'CnrdCnrd');
                a(Nw+ii,Nw+jj) = -0.5*Clah*chi05(ii)*dot(v,k_inf);
            end
        end

        % Canard-Wing influence
        for jj = 1:Nw
            v = computeHSV(ii,jj,P_w(jj,:),P_w(jj+1,:),P_h_mid(ii,:),ur,'CnrdWing');
            a(Nw+ii,jj) = -0.5*Clah*chi05(ii)*dot(v,k_inf);
        end
    end

    gamma   = a\b;
    gamma_w = gamma(1:Nw);
    gamma_h = gamma(Nw+1:end);

end