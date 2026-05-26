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