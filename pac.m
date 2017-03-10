global N m g Ix Iy Iz Jxz density S St r_CP r_CPt Cl0 Clalpha Cl0t Clalphat Cd0 K Cd0t Kt;

% Parameters
m = 0.800; % kg
density = 1.225; % kg/m^3
% I = [...]; 
S = 0.4; % m^2
St = 0.1; % m^2
r_CP = [-0.1;0;0]; % m
r_CPt = [-1;0;0]; % m

Cl0 = 0.01;
Clalpha = 2*pi;

Cl0t = 0;
Clalphat = 2*pi;

Cd0 = 0.05;
K = 0.1;

Cd0t = 0.02;
Kt = 0.9;

% Inertia
Ix = 0.2;
Iy = 0.2;
Iz = 0.2;
Jxz = 0.2;

% XFLR5
% Cl = @(alpha) 0.5*alpha + 0.01; % [rad]
% Cd = 0.1*Cl^2 + 0.01*Cl + 0.005;
% Cq = ...

% Contants
g = 9.80665; % m/s^2

% Lbw = [...];

% Solve ODE
tspan = [0 100]; % s

% global N;
N = 35;

y0 = zeros(N,1);
% y0(3) = -2;
% y0(34) = deg2rad(5);
% y0(7) = 8; % m/s

yp0 = zeros(N,1);

y0f = zeros(N,1);
% y0f(3) = 1;
% y0f(34) = 1;
% y0f(7) = 1;

yp0f = zeros(N,1);
% yp0f(2) = 1;
% yp0f(3) = 1;

% [y0n,yp0n] = decic(@odefun,0,y0,[],yp0,[]);

[tt,yy] = ode15i(@odefun,tspan,y0,yp0);

% Transform all equations to explicit and try ode45