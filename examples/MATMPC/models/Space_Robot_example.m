import casadi.*
clc
%% Import robot 
%--- URDF filename ---%

filename= 'R_RRRF.urdf';
%--- Create robot model ---%
[robot,robot_keys] = urdf2robot(filename);
n_q = robot.n_q;

%% Dimensions of state

nu_m = n_q; %No. of controls : torques of robot arm actuators

nx = 12 + 2*nu_m;  % Number of states - r0 + r0_dot + theta0 + omega0 + q + dq
%nu = 9;  % No. of controls
nz = 0;  % No. of algebraic states
ny = 18 + 3*nu_m; % Number of objective variables - r0 + r0_dot + theta0 + omega0 + q + qdot + omega0_dot + r0_ddot + q_ddot
% nyN = 27; 

%ny = 18; % Number of objective variables - p + p_dot + eta + omega + q + qdot 
nyN = 18; 
np = 0; % Number of parameters
nc = 0; % Number of general constraints
ncN = 0;

    % No. of controls: base
nu = nu_m;  % No. of controls

%% Constraints on state - number and indices
nbx = 11;
nbx_idx = [1:3, 7:8, 13:13+nu_m-1, 13+nu_m:13+2*nu_m-1];

%% Constraints on input - number and indices
nbu = nu;
nbu_idx = 1:nbu; % Index of control bounds 

%% Vectors for system equation
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

% Input variables definition

tau_q = SX.sym('controls', n_q, 1);

for i = 1:n_q
    tau_q(i) = controls(i);
end
% tau1 = controls(7);
% tau2 = controls(8);
% tau3 = controls(9);


%% --- States Parameters ---%
% base position and velocity (transl)
r0=[states(1); states(2); states(3)];
r0_dot=[states(4); states(5); states(6)];

% Euler angles and base frame angular vel 
phi = states(7);
theta = states(8);
psi = states(9);

theta0 = [phi; theta; psi];
R0 = Angles123_DCM(theta0);

omega0 = states(10:12);
u0 = [r0_dot ; omega0];
% robot arms joints positions
q       = states(13:13+n_q-1);
qdot    = states(13+n_q:end);
gen_vel = [u0; qdot];

%--- Kinematics ---%
%Kinematics
[RJ,RL,rJ,rL,e,g]=Kinematics(R0,r0,q,robot);
%Diferential Kinematics
[Bij,Bi0,P0,pm]=DiffKinematics(R0,r0,rL,e,g,robot);
%Velocities
[t0,tm]=Velocities(Bij,Bi0,P0,pm,u0,qdot,robot);
% %Jacobian of the last link
% [J0n, Jmn]=Jacob(rL(1:3,end),r0,rL,P0,pm,robot.n_links_joints,robot);
N    = NOC(r0,rL,P0,pm,robot);
Ndot = NOCdot(r0,t0,rL,tm,P0,pm,robot);
%--- Dynamics ---%
%Inertias in inertial frames
[I0,Im]=I_I(R0,RL,robot);

%Generalized Inertia matrix
H = GIM_NOC(N,I0,Im,robot); 
%Generalized Convective Inertia matrix
C = CIM_NOC(N,Ndot,t0,tm,I0,Im,robot);
% H = [H0, H0m; H0m.', Hm];
% C = [C0, C0m; Cm0  , Cm];

invH = simplify(inv(H)); %solve(M, SX.eye(M.size1())); % check inv in casadi in python api 
Qv = C*gen_vel;
real_u = [SX.zeros(6,1); tau_q];

epsilon_ddot = invH*(real_u-Qv);
%
epsilon_ddot = simplify(epsilon_ddot);

T_plate = [1 cos(phi)*tan(theta)  sin(phi)*tan(theta);
           0 cos(phi)            -sin(phi);
           0 sin(phi)/cos(theta)  cos(phi)/cos(theta)];     
       
r0_all  = r0_dot;
acc_all = epsilon_ddot(1:3); 
att_dot = T_plate*omega0;
domega  = epsilon_ddot(4:6);
dq      = qdot;
ddq     = epsilon_ddot(7:end); 

x_dot   = [
        r0_all;
        acc_all;
        att_dot;
        domega;
        dq;
        ddq;
        ];

xdot = SX.sym('xdot', nx, 1);
impl_f = xdot - x_dot;
%% Objectives and constraints
% Objective variables

h = [r0;r0_dot; theta0; omega0; q; dq ; acc_all; domega ; ddq]; 
% p + eta + p_dot + omega + q + q_dot +  p_ddot + omega_dot + q_ddot

%h = [x; y; z; phi; theta; psi; u; v; w; p; q; r; q1; q2; q3; dq]; 
% p + eta + p_dot + omega + q + q_dot    
hN = h(1:nyN);

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

Ts = 0.025;
Ts_st = 0.1; % shooting interval time