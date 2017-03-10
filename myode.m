function Yp = myode(t,Y)
    global N m g Ix Iy Iz Jxz density S St r_CP r_CPt Cl0 Clalpha Cl0t Clalphat Cd0 K Cd0t Kt;

    psi = Y(4);
    theta = Y(5);
    phi = Y(6);

    u = Y(7);
    v = Y(8);
    w = Y(9);

    p = Y(10);
    q = Y(11);
    r = Y(12);

    Fx = Y(13);
    Fy = Y(14);
    Fz = Y(15);

    Mx = Y(16);
    My = Y(17);
    Mz = Y(18);

    L = Y(19);
    D = Y(20);
    Q = Y(21);

    Cl = Y(22);
    Cd = Y(23);
    Cq = Y(24);

    Fxt = Y(25);
    Fyt = Y(26);
    Fzt = Y(27);

    Lt = Y(28);
    Dt = Y(29);
    Qt = Y(30);

    Clt = Y(31);
    Cdt = Y(32);
    Cqt = Y(33);

    alpha = Y(34);
    beta = Y(35);

    xp = Yp(1);
    yp = Yp(2);
    zp = Yp(3);
    psip = Yp(4);
    thetap = Yp(5);
    phip = Yp(6);
    up = Yp(7);
    vp = Yp(8);
    wp = Yp(9);
    pp = Yp(10);
    qp = Yp(11);
    rp = Yp(12);
%     Fxp = Yp(13);
%     Fyp = Yp(14);
%     Fzp = Yp(15);
%     Mxp = Yp(16);
%     Myp = Yp(17);
%     Mzp = Yp(18);
%     Lp = Yp(19);
%     Dp = Yp(20);
%     Qp = Yp(21);
%     Clp = Yp(22);
%     Cdp = Yp(23);
%     Cqp = Yp(24);
%     Fxtp = Yp(25);
%     Fytp = Yp(26);
%     Fztp = Yp(27);
%     Ltp = Yp(28);
%     Dtp = Yp(29);
%     Qtp = Yp(30);
%     Cltp = Yp(31);
%     Cdtp = Yp(32);
%     Cqtp = Yp(33);
%     alphap = Yp(34);
%     betap = Yp(35);
    
    Lbw = [ cos(alpha)*cos(beta)    -cos(alpha)*sin(beta)   -sin(alpha) ;
            sin(beta)               cos(beta)               0           ;
            sin(alpha)*cos(beta)    -sin(alpha)*sin(beta)   cos(alpha)  ];

    f = zeros(N,1);

    f(1) = -xp + (cos(theta)*cos(psi)) * u + (sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)) * v + (cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi)) * w;
    f(2) = -yp + (cos(theta)*sin(psi)) * u + (sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi)) * v + (cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi)) * w;
    f(3) = -zp + (-sin(theta) * u) + (sin(phi)*cos(theta)) * v + (cos(phi)*cos(theta)) * w;
    
    f(4)=-m*g*sin(theta)+Fx+Fxt-m*(up-r*v+q*w);
    f(5)=m*g*cos(theta)*sin(phi)+Fy+Fyt-m*(vp+r*u-p*w);
    f(6)=m*g*cos(theta)*cos(phi)+Fz+Fzt-m*(wp-q*u+p*v);
    f(7)=Mx-Ix*pp+Jxz*rp-(Iz-Iy)*q*r+Jxz*p*q;
    f(8)=My-Iy*qp+(Iz-Ix)*p*r-Jxz*(p^2-r^2);
    f(9)=Mz-Iz*rp-Iz*rp+Jxz*pp+(Ix-Iy)*p*q+Jxz*q*r;
    f(10)=p-phip+psip*sin(theta);
    f(11)=q-thetap*cos(phi)-psip*cos(theta)*sin(phi);
    f(12)=r+thetap*sin(phi)-psip*cos(theta)*cos(phi);
    
    f(13) = 1/2*density*(u^2+v^2+w^2)*S*Cl-L;
    f(14) = 1/2*density*(u^2+v^2+w^2)*S*Cd-D;
    f(15) = 1/2*density*(u^2+v^2+w^2)*S*Cq-Q;
    f(16) = 1/2*density*(u^2+v^2+w^2)*St*Clt-Lt;
    f(17) = 1/2*density*(u^2+v^2+w^2)*St*Cdt-Dt;
    f(18) = 1/2*density*(u^2+v^2+w^2)*St*Cqt-Qt;
    
    ff = -[Mx;My;Mz] + cross(r_CP,[Fx;Fy;Fz]) + cross(r_CPt,[Fxt;Fyt;Fzt]);
    f(19) = ff(1);
    f(20) = ff(2);
    f(21) = ff(3);
    
    ff = -[Fx;Fy;Fz] - Lbw * ([D;Q;L]);
    f(22) = ff(1);
    f(23) = ff(2);
    f(24) = ff(3);
    
    ff = -[Fxt;Fyt;Fzt] - Lbw * ([Dt;Qt;Lt]);
    f(25) = ff(1);
    f(26) = ff(2);
    f(27) = ff(3);
    
    f(28) = -Cl + Cl0 + Clalpha * alpha;
    f(29) = -Cd + Cd0 + K * Cl^2;
    f(30) = -Cq + Cd;

    f(31) = -Clt + Cl0t + Clalphat * alpha;
    f(32) = -Cdt + Cd0t + Kt * Clt^2;
    f(33) = -Cqt + Cdt;

%     ff = -[u;v;w] + Lbw * [norm([u;v;w]);0;0];
%     f(34) = ff(1);
%     f(35) = ff(2);
%     f(36) = ff(3);

end