%---------------------------------------------------------------%

% Model TiltHex considering the traslational part in world frame

%---------------------------------------------------------------%


%% Dimensions
nx = 18;  % Number of states - p + p_dot + eta + omega + forces 
nu = 6;  % No. of controls
nz=0;  % No. of algebraic states
ny = 18; % Number of objective variables - p + p_dot + eta + omega + omega_dot + p_ddot
nyN = 18; 
np = 0; % Number of parameters
nc = 0; % Number of general constraints
ncN = 0; 
nu_v = 0;       % No. of controls: servo_angle_rates
nu_lambda = 6;      % No. of controls: actuator_force_derivatives
nu = nu_lambda + nu_v;  % No. of controls

%% Constraints on state - number and indices
nbx = 11;
nbx_idx = [1:3, 7:8, 13:13+nu-1];

%% Constraints on input - number and indices
nbu = nu;
nbu_idx = 1:nbu; % Index of control bounds 
addpath('/home/ichelikh/Downloads/casadi_3.3.0')
import casadi.*

states   = SX.sym('states', nx, 1);
controls = SX.sym('controls', nu, 1);
alg      = SX.sym('alg',nz,1);
params   = SX.sym('paras', np, 1);
refs     = SX.sym('refs', ny, 1);
refN     = SX.sym('refs', nyN, 1);
Q        = SX.sym('Q', ny, 1);
QN       = SX.sym('QN', nyN, 1);
aux      = SX.sym('aux',ny,1);      % auxilary variable
auxN     = SX.sym('auxN',nyN,1);    % auxilary variable'

%% Variables definition
% State variables definition
x = states(1);
y = states(2);
z = states(3);
u = states(4);  
v = states(5);
w = states(6);
phi = states(7);
theta = states(8);
psi = states(9);
p = states(10);
q = states(11);
r = states(12);
f1 = states(13);
f2 = states(14);
f3 = states(15);
f4 = states(16);
f5 = states(17);
f6 = states(18);

% Input variables definition
df1 = controls(1);
df2 = controls(2);
df3 = controls(3);
df4 = controls(4);
df5 = controls(5);
df6 = controls(6);

%% Other variables definition

% Gravity [m/s^2]
g = 9.81;

% Mass of TiltHex [kg]
m = 2.771; % weight with the battery
indata.m = m;
% Inertia [kg*m^2]
Ix = 0.083;
Iy = 0.081;
Iz = 0.16;

% Propeller coefficients
c_f = 9.9e-4;
c_t = 1.8810e-05;
c = c_t / c_f; % f^c_tau in Davide's paper

% Allocation matrix
G = [  -0.2449    0.2673   -0.2832    0.2496    0.6017   -0.5546
        0.4993   -0.4947   -0.4691    0.4971    0.0212   -0.0038
        0.7157    0.6958    0.6572    0.6714    0.7051    0.7304
       -0.0031    0.1762    0.1769   -0.0046   -0.1827   -0.1746
       -0.2100   -0.0996    0.0945    0.2112    0.1111   -0.1077
        0.1976   -0.1857    0.1864   -0.1945    0.1894   -0.1923    ];

G1 = G(1:3, :);
G2 = G(4:6, :);

% Rotation Matrix Elements
r11 = cos(theta) * cos(psi);
r21 = cos(theta) * sin(psi);
r31 = -sin(theta);
r12 = sin(phi) * sin(theta) * cos(psi) - cos(phi) * sin(psi);
r22 = sin(phi) * sin(theta) * sin(psi) + cos(phi) * cos(psi);
r32 = sin(phi) * cos(theta);
r13 = cos(phi) * sin(theta) * cos(psi) + sin(phi) * sin(psi);
r23 = cos(phi) * sin(theta) * sin(psi) - sin(phi) * cos(psi);
r33 = cos(phi) * cos(theta);

R = [[r11 r21 r31].', [r12 r22 r32].', [r13 r23 r33].'];

% Body-wrench matrices
A1 = R * G1;
A2 = G2;

%% Dynamics
% Dynamic equations of the model
x_dot = [ 
          u;
          v;
          w;
          0       +         (A1(1,1)*f1 + A1(1,2)*f2 + A1(1,3)*f3 + A1(1,4)*f4 + A1(1,5)*f5 + A1(1,6)*f6)/m ;
          0       +         (A1(2,1)*f1 + A1(2,2)*f2 + A1(2,3)*f3 + A1(2,4)*f4 + A1(2,5)*f5 + A1(2,6)*f6)/m ;
          -g      +         (A1(3,1)*f1 + A1(3,2)*f2 + A1(3,3)*f3 + A1(3,4)*f4 + A1(3,5)*f5 + A1(3,6)*f6)/m ;
          p+r*cos(phi)*tan(theta)+q*sin(phi)*tan(theta);
          q*cos(phi) - r*sin(phi);
          r*cos(phi)/cos(theta) + q*sin(phi)/cos(theta);
          q*r*(Iy-Iz)/Ix +         (A2(1,1)*f1 + A2(1,2)*f2 + A2(1,3)*f3 + A2(1,4)*f4 + A2(1,5)*f5 + A2(1,6)*f6)/Ix ;
          r*p*(Iz-Ix)/Iy +         (A2(2,1)*f1 + A2(2,2)*f2 + A2(2,3)*f3 + A2(2,4)*f4 + A2(2,5)*f5 + A2(2,6)*f6)/Iy ;
          q*p*(Ix-Iy)/Iz +         (A2(3,1)*f1 + A2(3,2)*f2 + A2(3,3)*f3 + A2(3,4)*f4 + A2(3,5)*f5 + A2(3,6)*f6)/Iz ;
          df1;
          df2;
          df3;
          df4;
          df5;
          df6
    ];


xdot = SX.sym('xdot', nx, 1);
% Implicit f
impl_f = xdot - x_dot;
     
%% Objectives and constraints
% Objective variables
h = [x; y; z; phi; theta; psi; u; v; w; p; q; r; x_dot(4:6); x_dot(10:12)]; % p + eta + p_dot + omega + p_ddot + omega_dot
hN = h;

% Outer objectives
% Use diag(Q) or diag(QN) because Q and QN are vectors!
obji = 0.5 * (h - refs)' * diag(Q) * (h - refs);
objN = 0.5 * (hN - refN)' * diag(QN) * (hN - refN);

obji_GGN = 0.5*(aux-refs)'*(aux-refs);
objN_GGN = 0.5*(auxN-refN)'*(auxN-refN);
% General constraints
general_con = [];
general_con_N = [];

%% NMPC sampling time [s]

Ts = 0.004;
Ts_st = 0.1; % shooting interval time

%% Export
indata.model.g = g;
indata.model.m = m;
indata.model.G = G;
indata.model.J = diag([Ix Iy Iz]);
indata.model.c_f = c_f;
indata.model.c_t = c_t;
indata.model.nu_lambda = nu_lambda;
indata.model.nu_v = nu_v;
cd data
save('indata','indata');
cd ..