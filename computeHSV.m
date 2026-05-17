function v = computeHSV(ii,jj,PA,PB,P,ur,infl)

    rA = P - PA;
    rB = P - PB;

    nrA = norm(rA);
    nrB = norm(rB);
    urA = rA/nrA;
    urB = rB/nrB;
    
    if ii == jj && (strcmp(infl,'WingWing') || strcmp(infl,'CnrdCnrd'))
        vInfA = 1/(4*pi) * ((1)/(nrA + dot(ur,rA))) * cross(ur,urA);
        vInfB = 1/(4*pi) * ((1)/(nrB + dot(ur,rB))) * cross(ur,urB);
        v     = vInfA - vInfB;
    else
        vAB   = 1/(4*pi) * ((nrA + nrB)/(nrA*nrB*(nrA*nrB + dot(rA,rB)))) * cross(rA,rB);
        vInfA = 1/(4*pi) * ((1)/(nrA + dot(ur,rA))) * cross(ur,urA);
        vInfB = 1/(4*pi) * ((1)/(nrB + dot(ur,rB))) * cross(ur,urB);
        v     = vInfA + vAB - vInfB;
    end

end