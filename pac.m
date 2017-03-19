clear variables;
close all;

global T TT ZZ Z Alpha Lift Vel m g Ix Iy Iz Jxz density S c Cl0 Clalpha Cd0 K Cmx0 Cmy0 Cmz0 Cmxalpha Cmyalpha Cmzalpha Cmxbeta Cmybeta Cmzbeta;

global oT oY;

% Parameters
m = 0.800; % kg
density = 1.225; % kg/m^3
% I = [...]; 

S = 0.85; % m^2
c = 0.25; % m
% St = 0.1; % m^2
% r_CP = [-0.1;0;0]; % m
% r_CPt = [-1;0;0]; % m

Cl0 = 0.5;
Clalpha = 4.6;

% Cl0t = 0;
% Clalphat = 2*pi;

Cd0 = 0.015;
K = 1/(pi*S*0.9); % 0.03;

% Cd0t = 0.02;
% Kt = 0.9;

% Inertia
Ix = 0.45;
Iy = 0.39;
Iz = 0.84;
Jxz = 0.02;

% Moments
Cmx0 = 0;
Cmy0 = 0;
Cmz0 = 0;

Cmxalpha = 0; % -0.1;
Cmyalpha = -1;
Cmzalpha = 0; % 0.01;

Cmxbeta = 0; % 0.01;
Cmybeta = 0; % 0.1;
Cmzbeta = 0; % 0.05;

% XFLR5
% Cl = @(alpha) 0.5*alpha + 0.01; % [rad]
% Cd = 0.1*Cl^2 + 0.01*Cl + 0.005;
% Cq = ...

Alpha = 0;
Lift = 0;
Vel = 0;
T = 0;
TT = 0;
Z = 0;
ZZ = 0;

% Contants
g = 9.80665; % m/s^2

% Lbw = [...];

% Solve ODE
N = 12;
tspan = [0 100]; % s
y0 = zeros(N,1);
y0(3) = -2; % m
y0(7) = 8; % m/s

% Options
oT = 0;
oY = y0;
options = odeset('OutputFcn',@odeOutput); % ,'RelTol',1e-3,'AbsTol',1e-3);

[tt, yy] = ode45(@odeSystem, tspan, y0, options);

plot(yy(:,1),-yy(:,3));

% Number of variables/equations
% N = 35;
% 
% y0 = zeros(N,1);
% y0(3) = -2; % m
% y0(34) = deg2rad(5); % [rad]
% y0(7) = 8; % m/s
%     
% yp0 = zeros(N,1);
% 
% y0f = zeros(N,1);
% % y0f(1) = 1;
% % y0f(2) = 1;
% % y0f(3) = 1;
% % y0f(34) = 1;
% % y0f(7) = 1;
% 
% yp0f = zeros(N,1);
% yp0f(2) = 1;
% yp0f(3) = 1;
% 
% [y0n,yp0n] = decic(@odefun,0,y0,[],yp0,[]);
% 
% [tt,yy] = ode15i(@odefun,tspan,y0n,yp0n);

%%%

% Transform all equations to explicit and try ode45?

% NEW IDEA
% Solve algebraic equations first with values from Y
% and then subs the values found -> explicit ODEs
% ode45
