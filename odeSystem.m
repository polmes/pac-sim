function Yp = odeSystem(t,Y)
    global Z ZZ T TT Alpha Vel Lift m g Ix Iy Iz Jxz density S c Cl0 Clalpha Cd0 K Cmx0 Cmy0 Cmz0 Cmxalpha Cmyalpha Cmzalpha Cmxbeta Cmybeta Cmzbeta;
    
    global oT oY;
    
    T = [T; t];
    TT = [TT oT(end)];
    
    Z = [Z; Y(3)];
    ZZ = [ZZ; oY(3,end)];
    
%     x = Y(1);
%     y = Y(2);
%     z = Y(3);

    psi = Y(4);
    theta = Y(5);
    phi = Y(6);

    u = Y(7);
    v = Y(8);
    w = Y(9);

    p = Y(10);
    q = Y(11);
    r = Y(12);
    
    V = norm([u;v;w]);
    Vel = [Vel; V];
    
    beta = 0; % asin(v / V);
    alpha = atan2(w,u); % acos(u / (cos(beta) * V)); % should be negative?
    Alpha = [Alpha; alpha];
    
    Cl = Cl0 + Clalpha * alpha;
    Cd = Cd0 + K * Cl^2;
    Cq = 0; % Cd;
    
    L = 1/2*density*V^2*S*Cl;
    D = 1/2*density*V^2*S*Cd;
    Q = 1/2*density*V^2*S*Cq;
    Lift = [Lift; L];
    
    Lbw = [ cos(alpha)*cos(beta)    -cos(alpha)*sin(beta)   -sin(alpha) ;
            sin(beta)                cos(beta)               0          ;
            sin(alpha)*cos(beta)    -sin(alpha)*sin(beta)    cos(alpha) ];
    
    ff = Lbw * (-[D;Q;L]);
    Fx = ff(1);
    Fy = 0; % ff(2);
    Fz = ff(3);
    
    Cmx = 0; % Cmx0 + Cmxalpha * alpha + Cmxbeta * beta;
    Cmy = Cmy0 + Cmyalpha * alpha; %  + Cmybeta * beta;
    Cmz = 0; % Cmz0 + Cmzalpha * alpha + Cmzbeta * beta;
    
    Mx = 1/2*density*V^2*S*c*Cmx;
    My = 1/2*density*V^2*S*c*Cmy;
    Mz = 1/2*density*V^2*S*c*Cmz;
    
    xp = (cos(theta)*cos(psi)) * u + (sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)) * v + (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi)) * w;
    yp = (cos(theta)*sin(psi)) * u + (sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi)) * v + (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi)) * w;
    zp = (-sin(theta) * u) + (sin(phi)*cos(theta)) * v + (cos(phi)*cos(theta)) * w;
    
    psip = (q*sin(phi) + r*cos(phi)) * sec(theta);
    thetap = q*cos(phi) - r*sin(phi);
    phip = p + (q*sin(phi) + r*cos(phi)) * tan(theta);
    
    up = (Fx - m*g*sin(theta) - m*(-r*v + q*w)) / m;
    vp = (Fy + m*g*cos(theta)*sin(phi) - m*(r*u - p*w)) / m;
    wp = (Fz + m*g*cos(theta)*cos(phi) - m*(-q*u + p*v)) / m;
    
    pp = (Iz*Mx + Jxz*Mz + Jxz * ((Ix - Iy)*p*q - Jxz*q*r) - Iz * ((Iz - Iy)*q*r - Jxz*p*q)) / (Ix*Iz - Jxz^2);
    qp = (My + (Iz - Ix)*p*r - Jxz * (p^2-r^2)) / Iy;
    rp = (Jxz*Mx + Ix*Mz + Ix * ((Ix - Iy)*p*q - Jxz*q*r) - Jxz * ((Iz - Iy)*q*r - Jxz*p*q)) / (Ix*Iz - Jxz^2);
    
    Yp = [xp; yp; zp; psip; thetap; phip; up; vp; wp; pp; qp; rp];
end
