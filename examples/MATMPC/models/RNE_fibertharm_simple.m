%function RNE()
clc
%clear all
close all
% 
% addpath('/home/ichelikh/Downloads/casadi_3.3.0')
% addpath('/home/ichelikh/ms2021-idriss_chelikh/Software/robotics')
addpath('C:\Users\idris\Desktop\LAAS_AM\robotics')
addpath('C:\Users\idris\casadi-matlabR2014b-v3.3.0')
import casadi.*

%% Dimensions

nu_m = 3; %No. of controls : torques of robot arm actuators

nx = 12 + 2*nu_m;  % Number of states - p + p_dot + eta + omega + q + dq
%nu = 9;  % No. of controls
nz = 0;  % No. of algebraic states
ny = 18 + 3*nu_m; % Number of objective variables - p + p_dot + eta + omega + q + qdot + omega_dot + p_ddot + q_ddot
% nyN = 27; 

%ny = 18; % Number of objective variables - p + p_dot + eta + omega + q + qdot 
nyN = 18; 
np = 0; % Number of parameters
nc = 0; % Number of general constraints
ncN = 0;

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
q1 = states(13);
q2 = states(14);
q3 = states(15);
dq1 = states(16);
dq2 = states(17);
dq3 = states(18);
% f1 = states(19);
% f2 = states(20);
% f3 = states(21);
% f4 = states(22);
% f5 = states(23);
% f6 = states(24);

% Input variables definition

% f1 = controls(1);
% f2 = controls(2);
% f3 = controls(3);
% f4 = controls(4);
% f5 = controls(5);
% f6 = controls(6);
tau1 = controls(7);
tau2 = controls(8);
tau3 = controls(9);

%% ------------------------------MODELLING--------------------------------
%
% ------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Define the DH parameters

N_DOFS = nu_m;
dh.theta = SX.sym('theta'  , N_DOFS, 1);
dh.alpha = [pi/2 0 0];
dh.offset = [0 0 0];
dh.d = [0 0 0];
dh.a = [0 0.3 0.3];
dh.type = ['r' 'r' 'r'];



% -------------------------------------------------------------------------
% Rigid body paramaters: inertia, mass, and center of mass
I1 = diag([0.0005  0.000011 0.0004]);
I2 = diag([0.00035 0.0018   0.0016]);
I3 = diag([0.423   1.414    1.185])*1e-4;
rb.I =  repmat(eye(3), 1, 1, N_DOFS);

rb.I(:,:,1) = I1;
rb.I(:,:,2) = I2;
rb.I(:,:,3) = I3;

%rb.I =  repmat(zeros(3,3), 1, 1, N_DOFS);

%rb.m = [1 1 1];
rb.m = [0.059 0.256 0.093];
rb.mdrone = 2.122+0.512;
rb.I_drone = diag([0.0115,  0.0114, 0.0194]);
% In standard DH, COM is mesured respect to its end-effector (using local
% frame). When COM is [0 0 0]', it means the COM is located exactly at the
% end-effector. Therefore, COM usually has negative values, which means the
% COM is behind the end-effector
rb.r = [0 0 0; -0.15 0 0; -0.15 0 0]';
rb.rc0 = [0;0;-0.051];
rb.rc1 = [0;0;0];
 
%---------------------------------------------------------------------
% Here we go!

Gdrone = [ -0.2449    0.2673   -0.2832    0.2496    0.6017   -0.5546
            0.4993   -0.4947   -0.4691    0.4971    0.0212   -0.0038
            0.7157    0.6958    0.6572    0.6714    0.7051    0.7304
           -0.0031    0.1762    0.1769   -0.0046   -0.1827   -0.1746
           -0.2100   -0.0996    0.0945    0.2112    0.1111   -0.1077
            0.1976   -0.1857    0.1864   -0.1945    0.1894   -0.1923 ];

epsilon = [states(1:3); states(7:9); states(13:13+nu_m-1)]; % eps = [xb, phib, q]

depsilon = [states(4:6); states(10:12); states(13+nu_m : 13+2*nu_m-1)]; % deps = [vb, dphib, dq]

forces = controls(1:6);

phib = epsilon(4:6)';
wRb = att2RotMat('xyz', phib);

Gdrone1 = [wRb*Gdrone(1:3,:);
          Gdrone(4:6,:)];
           

Gforces = Gdrone1*forces;
Big_Eye = eye(nu);

tic
G_C_Mat = invdyn(dh, rb, epsilon, depsilon, zeros(nu,1),[0; 0; -9.8]);

M = SX.zeros(9,9);

for i = 1:nu
    M(:,i) = invdyn(dh, rb, epsilon, zeros(nu,1), Big_Eye(i,:)');
end 

M = simplify(M);

invM = simplify(inv(M)); %solve(M, SX.eye(M.size1())); % check inv in casadi in python api 

real_u = [Gforces;tau1;tau2;tau3];

epsilon_ddot = invM*(real_u-G_C_Mat);
%
epsilon_ddot = simplify(epsilon_ddot);

dh.theta = zeros(3,1);

T_plate = [1 cos(phi)*tan(theta)  sin(phi)*tan(theta);
           0 cos(phi)            -sin(phi);
           0 sin(phi)/cos(theta)  cos(phi)/cos(theta)];     
       
v_all   = [u;v;w];
acc_all = epsilon_ddot(1:3); 
att_dot = T_plate*[p;q;r];
domega  = epsilon_ddot(4:6);
dq      = [dq1; dq2; dq3];
ddq     = epsilon_ddot(7:9); 


x_dot   = [
        v_all;
        acc_all;
        att_dot;
        domega;
        dq;
        ddq;
        ];
    
toc
xdot = SX.sym('xdot', nx, 1);
impl_f = xdot - x_dot;
%% Objectives and constraints
% Objective variables

 h = [x; y; z; phi; theta; psi; u; v; w; p; q; r; q1; q2; q3; dq ; acc_all; domega ; ddq]; 
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

%% Other physical parameters for simulation
c_f = 9.9e-4;
c_t = 1.8810e-05;
c = c_t / c_f; % f^c_tau in Davide's paper


%% Export
%indata.model.dh = dh;
indata.model.rb = rb;
indata.model.G = Gdrone;

indata.model.c_f = c_f;
indata.model.c_t = c_t;
indata.model.nu_lambda = nu_lambda;
indata.model.nu_v = nu_v;

cd data
save('indata','indata');
cd ..
%% -----------------FUNCTIONS USED FOR MODEL COMPUTATION-------------------

% -------------------------------------------------------------------------
function Qtot = invdyn(dh, rb, epsilon, depsilon, ddepsilon,  grav)
% Inverse dynamic with recursive Newton-Euler following De luca slides
import casadi.*

if nargin < 6
    grav = [0;0;0];
end

z0 = [0; 0; 1];

bwb    = depsilon(4:6);  %omega of base com in base frame
wdotb = ddepsilon(4:6); %omega dot of base com in base frame
a0b   = ddepsilon(1:3); %p ddot of base com in world frame

q     = epsilon(7:end); 
qdot  = depsilon(7:end);
qddot = ddepsilon(7:end);


N_DOFS = length(q);
% init all to sx zeros
w     = SX.sym('w'   , 3,N_DOFS);
wdot  = SX.sym('wdot', 3,N_DOFS);
vdot  = SX.zeros(3,N_DOFS);
n     = SX.sym('n', 3,N_DOFS);
z     = SX.zeros(3,N_DOFS);
% F     = SX.sym('F', 3,1);
% N     = SX.sym('N', 3,1);

f     = SX.sym('f', 3,N_DOFS);

Q     = SX.sym('Q', N_DOFS,1);

%Rb = eye(3);

R0b = [0 1 0;
       0 0 1;
       1 0 0];

%init with respect of world frame

phib = epsilon(4:6)'; %yaw-pitch-roll Euler angles (moving)
wRb  = att2RotMat('xyz', phib);
p0b  = [0; 0; -0.05];  % offset between link 0 and com in BASE frame

% omega and omega dot of base in world frame
wwb    = wRb*bwb;
wdotwb = wRb*wdotb;
%acc resultant of 1st link frame 0 with respect of worldframe
w_a0 = a0b + cross(wdotwb , (wRb*p0b)) + cross(wwb, cross(wwb, wRb*p0b)) - grav; % premulitply with rotations
% you need to rotate in frame 0 the gravity (which is expressed in world frame) 
a0 = R0b*wRb'*w_a0;
w0 = R0b*bwb; 
wdot0 = R0b*wdotb;
    % ---------------------------------------------------------------------
    % Forward recursion 
    for i = 1 : N_DOFS
        T = calc_transformation(i-1, i, dh, q);
        R = T(1:3, 1:3);
        %R0 = R0*R;
        p = [dh.a(i); dh.d(i)*sin(dh.alpha(i)); dh.d(i)*cos(dh.alpha(i))]; % opposite of formula "La matrice de Denavit-Hartenberg modifiée : [..]" at https://fr.wikipedia.org/wiki/Denavit-Hartenberg
        %z = R*z0;
        if i > 1
            
            z(:,i)    =  z0; %R*z(:,i-1);
            w(:,i)    =  R'*(w(:,i-1) + z(:,i-1).*qdot(i)); % substitute z(:,i-1) with z(:,i) <============
             
            wdot(:,i) =  R'*(wdot(:,i-1) +  z(:,i-1).*qddot(i) + ... % 2x substitute z(:,i-1) with z(:,i)<============
                         cross(w(:,i-1), z(:,i-1).*qdot(i)));
                    
            vdot(i,:) =  R'*vdot(:,i-1) + cross(wdot(:,i), p) + ...
                         cross(w(:,i), cross(w(:,i),p));
        else
            z(:,i)    =  z0;  %R*z0;
            w(:,i)    =  R'*(w0 + z0.*qdot(i)); % b_w_b = 0_w_0 = w0 not true because of initial pi offset on q1
            
            wdot(:,i) =  R'*((wdot0 + z0.*qddot(i))+ ... 
                         cross(w0, z0.*qdot(i)));
            vdot(:,i) =  R'*(a0) + cross(wdot(:,i), p) + ...
                         cross(w(:,i), cross(w(:,i),p));
        end
    end
    
    z0 = [0; 0; 1];
    
    % Dynamic simulation
    % Backward recursion
    for i = N_DOFS:-1:1
        p = [dh.a(i); dh.d(i)*sin(dh.alpha(i)); dh.d(i)*cos(dh.alpha(i))];
        
        vcdot = vdot(:,i) + cross(wdot(:,i),rb.r(:,i)) + ...
                cross(w(:,i),cross(w(:,i),rb.r(:,i)));
        
        F = rb.m(i)*vcdot;
        N = rb.I(:,:,i)*wdot(:,i)+cross(w(:,i),rb.I(:,:,i)*w(:,1));
        
        if i < N_DOFS
            T = calc_transformation(i, i+1, dh, q);
            R = T(1:3, 1:3);
            
            n(:,i) = R*(n(:,i+1) + cross(R'*p, f(:,i+1))) + ...
                     cross(rb.r(:,i)+p,F) + N;
                 
            f(:,i) = R*f(:,i+1) + F;
        else
            n(:,i) = cross(rb.r(:,i)+p,F) + N;
            f(:,i) = F;
        end
        
        T = calc_transformation(i-1, i, dh, q);
        R = T(1:3, 1:3);
      % z0 = [0;0;1];
        if dh.type(i) == 't'
            Q(i) = f(:,i)'*R'*z(:,i) ;
        elseif dh.type(i) == 'r'
            Q(i) = n(:,i)'*R'*z(:,i) ;
        end
    end
    
    wrenchB = wrench_base(f(:,1), n(:,1), R0b, p0b, wwb, rb.mdrone, rb.I_drone, ddepsilon, grav);
    
    Qtot = [wrenchB;Q];
    
end

function [wrench] = wrench_base(f0, n0, Rb0, p0b, bwb, m, J, ddepsilon, grav)

%   RotWM = [Rw0,      zeros(3,3);
%            zeros(3,3), eye(3)];
  RotWM = eye(6);

  Mass  = [ m*eye(3),  zeros(3,3); 
            zeros(3,3),        J];
  
  RotM0 = [[eye(3)         , zeros(3,3)];
           [skewMat(p0b)   , Rb0       ]];
  f0 = Rb0'*f0;
  n0 = Rb0'*n0;
  wrench0 = [f0; n0];
  
  gyro = -cross(bwb, J*bwb);

  Grav_and_Gyro = [m*grav ; gyro];
  
  acc = ddepsilon(1:6);

  wrench =  RotWM*(Mass*acc - Grav_and_Gyro - RotM0*wrench0);
end
% -------------------------------------------------------------------------
function  T = calc_transformation(from, to, dh, qc)
% Transformation from one joint to another joint
% 0<=from<N_DOFS
% 0<to<=N_DOFS

T = eye(4);
N_DOFS = length(qc);

% Sanity check
if (from >= N_DOFS) || (from < 0) || (to <= 0) || (to >  N_DOFS)
    return;
end

for i = from+1 : to
    if     dh.type(i) == 'r'
           dh.theta(i) = qc(i);
    elseif dh.type(i) == 'p'
           dh.d(i) = qc(i);
    end

    ct = cos(dh.theta(i));
    st = sin(dh.theta(i));
    ca = cos(dh.alpha(i));
    sa = sin(dh.alpha(i));
    
    T = T * [ ct    -st*ca   st*sa     dh.a(i)*ct ; ...
              st    ct*ca    -ct*sa    dh.a(i)*st ; ...
              0     sa       ca        dh.d(i)    ; ...
              0     0        0         1          ];
  
end
end
