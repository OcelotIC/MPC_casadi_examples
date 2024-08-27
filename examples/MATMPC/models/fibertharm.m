%---------------------------------------------------------------%

% Model TiltHex considering the traslational part in world frame

%---------------------------------------------------------------%


%% Dimensions
nx = 24;  % Number of states - p + p_dot + eta + omega + forces 
%nu = 9;  % No. of controls
nz=0;  % No. of algebraic states
ny = 18; % Number of objective variables - p + p_dot + eta + omega + omega_dot + p_ddot
nyN = 18; 
np = 0; % Number of parameters
nc = 0; % Number of general constraints
ncN = 0;
nu_m = 3; %No. of controls : torques of robot arm actuators
nu_v = 0;       % No. of controls: servo_angle_rates
nu_lambda = 6;      % No. of controls: actuator_force_derivatives
nu = nu_lambda + nu_v+nu_m;  % No. of controls

%% Constraints on state - number and indices
nbx = 17;
nbx_idx = [1:3, 7:8, 13:15, 16:18, 19:19+nu_lambda-1];

%% Constraints on input - number and indices
nbu = nu;
nbu_idx = 1:nbu; % Index of control bounds 
addpath('/home/ichelikh/Downloads/casadi_3.3.0')
addpath('/home/ichelikh/ms2021-idriss_chelikh/Software/robotics')
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
q1 = states(13);
q2 = states(14);
q3 = states(15);
dq1 = states(16);
dq2 = states(17);
dq3 = states(18);
f1 = states(19);
f2 = states(20);
f3 = states(21);
f4 = states(22);
f5 = states(23);
f6 = states(24);

% Input variables definition

df1 = controls(1);
df2 = controls(2);
df3 = controls(3);
df4 = controls(4);
df5 = controls(5);
df6 = controls(6);


% f1   = controls(1);
% f2   = controls(2);
% f3   = controls(3);
% f4   = controls(4);
% f5   = controls(5);
% f6   = controls(6);
tau1 = controls(7);
tau2 = controls(8);
tau3 = controls(9);


%% Other variables definition

robot.n = 3; % Number of links/joints (*assumption*: serial manipulator)
             % Only n firsts elements of the following vectors might be
             % used: thus define next vectors accordingly!
robot.joint_types = ['r' 'r' 'r']; % type of joint j. Valid values are: 'r': rotational, 'p':prismatic

%% ~~~~~~~~~~~~~~~~~~~~~~~~ DH PARAMETERS ~~~~~~~~~~~~~~~~~~~~~~~~
robot.dh.a = [0 0.3 0.3];
robot.dh.d = [0 0 0];
robot.dh.alpha = [pi/2 0 0];
robot.dh.theta = [q1 q2 q3];

% Values of the DH parameters when the robot is in its home (rest) pose
robot.dh.theta0 = [pi 0 0];

% Derivatives of the DH parameters
robot.dh.dtheta = [dq1 dq2 dq3];

%% ~~~~~~~~~~~~~~~~~~~~~~~~ KINEMATICS ~~~~~~~~~~~~~~~~~~~~~~~~
% pl: position of the c.o.m of Link i wrt Link Frame i (knwon at priori, depend upon robot's mechanical design)
robot.values_pl = [[0; 0; 0] [-0.3/2; 0; 0] [-0.3/2; 0; 0]]; 

% Position and orientation of Frame B wrt World Frame W
robot.base.pwb = [x; y; z];
robot.base.phib = [psi; theta; phi]; %yaw-pitch-roll Euler angles (moving)
robot.base.Rwb = att2RotMat('zyx', robot.base.phib);

% Linear and Angular velocities of Body wrt World Frame W
robot.base.dpwb = [u; v; w];
robot.base.dphib = [p; q; r];

% Position and orientation of Manipulator's Frame 0 wrt Body Frame B
%robot.base.pb0 = [xb0 yb0 zb0]';
robot.base.pb0 = [0 0 -0.1]';
robot.base.Rb0 = att2RotMat('zyx', [-pi/6 -pi/2 -pi/2]);


robot.base.Twb = [ cos(theta)*cos(psi)  -sin(psi)  0;
                   cos(theta)*sin(psi)   cos(psi)  0;
                        -sin(theta)        0       1 ];
robot.base.Qwb = simplify((robot.base.Rwb') * robot.base.Twb);

% Generalised joint vector
robot.epsilon = [ robot.base.pwb;
                  robot.base.phib;
                  robot.dh.theta' ];
                
robot.depsilon = [ robot.base.dpwb;
                   robot.base.dphib;
                   robot.dh.dtheta' ];

%% Compute Homogenous Transformation Matrices From Link Frame i-1 to Link Frame i
A = HTML2L(robot.dh.a, robot.dh.d, robot.dh.alpha, robot.dh.theta);
% A contains the homogeneous transformation matrices from link i-1 to link
% i.
% In the case of the 3dof arm, i is in {1, 2, 3}, where i = 3 is the
% end-effector "link", and Frame 0 denotes the "base" of the manipulator.
                 
%% Compute Links quantities
disp('');
disp('Preliminary kinematic quantities:');

m = 2.771; % weight with the battery
indata.m = m;

g = 9.81;

% Inertia [kg*m^2]
Ix = 0.083;
Iy = 0.081;
Iz = 0.16;

Hb = diag([Ix, Iy, Iz]); % inertial matrix of drone 
mb = indata.m; %kg mass of drone

H1 = diag([0.03, 0.03, 0.2]);
H2 = diag([0.03, 0.03, 0.2]);
H3 = diag([0.03, 0.03, 0.2]);

H_arm = zeros(3,3,robot.n);%inertia matrixes of arms bodies

H_arm(:,:,1)= H1;
H_arm(:,:,2)= H2;
H_arm(:,:,3)= H3;


m_arm = [0.5 0.5 0.3]; %masses of arm bodies
for i = 0 : robot.n
  
  if i == 0 % ------------ { base }

    robot.base.H = Hb;
    robot.base.m = mb;
    
    % Z-axis coordinates of Base Frame 0 wrt Body Frame B
    robot.base.zb0 = robot.base.Rb0 * [0; 0; 1];
  
  else % ------------ { link }
    fprintf('* for link %d\n', i);
    % Inertia (expressed wrt Link Frame i) and mass of Link i

    robot.links{i}.H = H_arm(:,:,i);
    robot.links{i}.m = m_arm(i);

    % Get Homogeneous Transformation Matrix from Link Frame i-1 to Link
    % Frame i
    A_i = A(1+(i-1)*4:4+(i-1)*4,:);
    % Get position and orientation of Link Frame i wrt Link Frame i-1 (DH
    % Parameters)
    robot.links{i}.R = A_i(1:3, 1:3);
    robot.links{i}.p = A_i(1:3, 4);
    
    % Position of the c.o.m of Link i wrt Link Frame i
      robot.links{i}.pl = robot.values_pl(:,i);
    
    % Position and orientation of Link Frame i wrt Link Base Frame 0
    if i == 1 % ------------ { first link }
      robot.links{i}.p0 = robot.links{i}.p;
      robot.links{i}.R0 = robot.links{i}.R;

    else % ------------ { i-th link }
      robot.links{i}.p0 = robot.links{i-1}.p0 + robot.links{i-1}.R0 * robot.links{i}.p;
      robot.links{i}.R0 = robot.links{i-1}.R0 * robot.links{i}.R;
      
    end
    
    % Position and orientation of Link Frame i wrt Body Frame B
    robot.links{i}.pb = robot.base.pb0 + robot.base.Rb0 * robot.links{i}.p0;
  
    robot.links{i}.Rb = robot.base.Rb0 * robot.links{i}.R0;

    % Z-axis coordinates of Link Frame i wrt Base Frame 0
    %   [0; 0; 1] = zi 
    %   zi: z-axis coordinates of Link Frame i wrt Link Frame i
    robot.links{i}.z0 = robot.links{i}.R0 * [0; 0; 1]; 
    
    % Z-axis coordinates of Link Frame i wrt Body Frame B
    robot.links{i}.zb = robot.base.Rb0 * robot.links{i}.z0;

    % Position of the c.o.m of Link i wrt Body Frame B
    robot.links{i}.pbl = robot.links{i}.pb + robot.links{i}.Rb * robot.links{i}.pl;
  end
         
end

%% Kinematics
disp('');
disp('Kinematics:')

for i = 1:robot.n % cycle over links
  fprintf('* Jacobian of link %d\n', i);
  
  %Ji = sym('j%d%d', [6 robot.n], 'real'); 
  Ji = SX.sym('J',6, robot.n);
  for j = 1:robot.n % ------------ { cycle over joints (columns of each Link i's Jacobian) }

    if j <= i % ------------ { Effects of Links before i or i itself }

      if j == 1 % ------------ { first joint }
        switch lower(robot.joint_types(j))
          case 'p'
            Ji(:,j) = [ robot.base.zb0;
                        zeros(3,1)
                      ];
          case 'r'
            Ji(:,j) = [ skewMat(robot.base.zb0)*(robot.links{i}.pbl - robot.base.pb0); %issue w/ cross product double cross SX error
                        robot.base.zb0                      
                      ];
          otherwise
            fprintf('Joint %d\n', j); error('Invalid joint type');
        end

      else % ------------ { j-th joint }
        switch lower(robot.joint_types(j))
          case 'p'
            Ji(:,j) = [ robot.links{j-1}.zb;
                        zeros(3,1)
                      ];
          case 'r'
            Ji(:,j) = [ skewMat(robot.links{j-1}.zb) * (robot.links{i}.pbl - robot.links{j-1}.pb);
                        robot.links{j-1}.zb                     
                       ];
          otherwise
            fprintf('Joint %d\n', j); error('Invalid joint type');
        end
      end

    else % ------------ { Effects of Links after i }
      Ji(:,j) = zeros(6, 1);
    end
  end
  
  Si = skewMat(robot.base.Rwb * robot.links{i}.pbl);
  Hbli = robot.links{i}.Rb * robot.links{i}.H * (robot.links{i}.Rb');
  
  robot.links{i}.J    = Ji;
  robot.links{i}.S    = Si;
  robot.links{i}.Hbl  = Hbli;
end 

%% Dynamics
disp('');
disp('Dynamics:')

e3 = [0; 0; 1];
epsilon = robot.epsilon;
depsilon = robot.depsilon;
nepsilon = length(robot.epsilon);

% Get variables related to the Base
mb    = robot.base.m;
Hb    = robot.base.H;
pb    = robot.base.pwb;
Rb    = robot.base.Rwb;
Tb    = robot.base.Twb;
Qm     = robot.base.Qwb;

% Initialize submatrices of Generalised Inertia Matrix B
B11 = mb*eye(3);

B22 = (Qm') * Hb * Qm;

B33 = zeros(3,3);

B12 = zeros(3,3);

B13 = zeros(3,3);

B23 = zeros(3,3);

% Initialize Potential Energy U
U = mb * g * (e3') * pb;


disp('* Generalised Inertia Matrix B')

for i = 1:robot.n
  
  fprintf('** Link %d\n', i);
  
  % Get variables related to Link i
  Jpi   = robot.links{i}.J(1:3,:);
  Joi   = robot.links{i}.J(4:6,:);
  mi    = robot.links{i}.m;
  Si    = robot.links{i}.S;
  Hbli  = robot.links{i}.Hbl;
  pbli  = robot.links{i}.pbl;
  
  B11 = B11 + mi;
  
  B22 = B22 ...
            + mi * (Tb)' * (Si') * Si * Tb ...
            + (Qm') * Hbli * Qm;
  
  B33 = B33 ...
            + mi * (Jpi') * Jpi ...
            + (Joi') * Hbli * Joi;
        
  B12 = B12 ...
            - mi * (Si) * Tb;
        
  B13 = B13 ...
            + mi * Rb * Jpi;
  
  B23 = B23 ...
            + (Qm') * Hbli * Joi ...
            - mi * (Tb') * (Si') * Rb * Jpi;
          
  U = U ...
        + mi * g * (e3') * (pb + Rb * pbli);
          
end

% Finalize submatrices
B11 = B11 * eye(3);


% Exploit symmetry of Generalised Inertia Matrix
B21 = B12';
B31 = B13';
B32 = B23';

% Compose Generalised Inertia Matrix B
B = [ B11 B12 B13;
      B21 B22 B23;
      B31 B32 B33 ];

disp('* Gravitational Vector G')
G = jacobian(U, epsilon);

%% compute C matrix Coriolis
C = SX.sym('C', nepsilon,nepsilon);
for i = 1:nepsilon
  
  epsiloni = epsilon(i);
  
  % Fill all columns at row i of Coriolis Matrix C
  % jacobian(B,epsiloni)
  Bk = [];
  for k =1:nepsilon
    Bk = [Bk jacobian(B(:, k),epsiloni)];
    
  end
    Ctemp = 0.5 * (jacobian(B(:,i),epsilon) + jacobian(B(:,i),epsilon).' - Bk);
    Ctemp = depsilon.' * Ctemp * depsilon;
    C(:,i) = Ctemp; % here you get C*epsilon, which is different from C !!
end
%% allocation matrix of drone base 

Gdrone = [-0.0139,  0.2777, -0.2991, -0.0109,  0.2729, -0.2964;
          0.3343, -0.1620, -0.1585,  0.3330, -0.1601, -0.1595;
          0.9486,  0.9360,  0.9494,  0.9391,  0.9514,  0.9389;
          0.0056,  0.3021,  0.3042,  0.0044, -0.2930, -0.2942;
         -0.3525, -0.1784,  0.1640,  0.3383,  0.1633, -0.1782;
          0.1247, -0.1661,  0.1234, -0.1681,  0.1227, -0.1668];


G1 = Gdrone(1:3, :);
G2 = Gdrone(4:6, :);

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

G1 = R * G1;
G2 = G2;
G_drone = [G1;
           G2];
Q_nc = [G_drone,zeros(6, 3) ;
        zeros(3, 6), eye(3)];

%diag(A_drone, eye(nu_m));
%% EoM inversion 
 invB = solve(B, SX.eye(B.size1()));
% 
forces = [f1 ;f2; f3; f4; f5; f6; tau1; tau2; tau3];
 ddeps = invB*(Q_nc*forces - C*depsilon - G');

v_all   = [u;v;w];
acc_all = ddeps(1:3); 
att_dot     = [p;q;r];
domega  = ddeps(4:6);
dq      = [dq1; dq2; dq3];
ddq     = ddeps(7:9); 
dforce  = [df1;df2;df3;df4;df5;df6];

x_dot   = [
        v_all;
        acc_all;
        att_dot;
        domega;
        dq;
        ddq;
        dforce
        ];
    
xdot = SX.sym('xdot', nx, 1);
% Implicit f
impl_f = xdot - x_dot;

%% Objectives and constraints
% Objective variables
h = [x; y; z; phi; theta; psi; u; v; w; p; q; r; q1; q2; q3; dq]; % p + eta + p_dot + omega + q + q_dot/ p_ddot + omega_dot not added
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

c_f = 9.9e-4;
c_t = 1.8810e-05;
c = c_t / c_f; % f^c_tau in Davide's paper
G = G';

%% Export
indata.model.g = g;
indata.model.m = m;
indata.model.G = Gdrone;
indata.model.J = diag([Ix Iy Iz]);
indata.model.c_f = c_f;
indata.model.c_t = c_t;
indata.model.nu_lambda = nu_lambda;
indata.model.nu_v = nu_v;

cd data
save('indata','indata');
cd ..