%function RNE()
clc
clear all
close all

addpath('/home/ichelikh/Downloads/casadi_3.3.0')
import casadi.*

% -------------------------------------------------------------------------
% Define the DH parameters
N_DOFS = 3;
dh.theta = SX.sym('theta'  , 3, 1);
dh.alpha = [0 0 0];
dh.offset = [0 0 0];
dh.d = [0 0 0];
dh.a = [0.3 0.3 0.3];
dh.type = ['r' 'r' 'r'];

% -------------------------------------------------------------------------
% Rigid body paramaters: inertia, mass, and cener of mass
rb.I =  repmat(eye(3), 1, 1, N_DOFS);
rb.m = [1 1 1];
% In standard DH, COM is mesured respect to its end-effector (using local
% frame). When COM is [0 0 0]', it means the COM is located exactly at the
% end-effector. Therefore, COM usually has negative values, which means the
% COM is behind the end-effector
rb.r = [-0.25 0 0; -0.25 0 0; -0.25 0 0]';

% -------------------------------------------------------------------------


%instead of trajectory, we use symbolics

q   = SX.sym('q'  , 3, 1);
dq  = SX.sym('dq' , 3, 1);
ddq = SX.sym('ddq', 3, 1);
tau = SX.sym('tau', 3, 1);
% -------------------------------------------------------------------------
% Kinematic simulation, optional, for visualization purpose!   
% End-effector position, form the base which is located at [0 0 0]'
EE = zeros(3, N_DOFS+1);
COM = zeros(3, N_DOFS);    


%% useless simulation of defined trqjectory, will be removed
% for k = 1 : length(time_span)  
%     for i = 1 : 1 : N_DOFS
%         T = calc_transformation(0, i, dh, qc(k,:));
%         EE(:,i+1) = T(1:3,4);
%         COM(:,i) = EE(:,i+1) + T(1:3,1:3) * rb.r(:,i);
%     end
%     
%     % Draw the robot
%     set(h1, 'XData', EE(1, :), 'YData', EE(2, :),'ZData', EE(3, :));
%     set(h2, 'XData', EE(1, :), 'YData', EE(2, :),'ZData', EE(3, :));
%     set(h3, 'XData', COM(1, :), 'YData', COM(2, :),'ZData', COM(3, :));
%     drawnow;
% end

%% -------------------------------------------------------------------------
% Here we go!
% Q = invdyn(dh, rb, qc, qcdot, qcddot, [0; -9.8; 0]);

w0 = SX.sym('w0'   , 3,1);
wdot0 = SX.sym('wdot0'   , 3,1);
vdot0 = SX.sym('vdot0'   , 3,1);
null = zeros(3,1);
tau = invdyn(dh, rb, q, dq, ddq, [0; 0;-9.8], tau, w0, wdot0, vdot0 );

Grav_term = invdyn(dh, rb, q, zeros(3,1), zeros(3,1), [0; 0; -9.8], tau, null, null, null)
% M1 = invdyn(dh, rb, q, zeros(3,1), [1; 0; 0], [0; 0; 0], tau, null, null, null);
% M2 = invdyn(dh, rb, q, zeros(3,1), [0; 1; 0], [0; 0; 0], tau, null, null, null);
% M3 = invdyn(dh, rb, q, zeros(3,1), [0; 0; 1], [0; 0; 0], tau, null, null, null);
% CandM = invdyn(dh, rb, q, dq, zeros(3,1), [0; 0; -9.8], tau, null, null, null);
% 
% M = [M1 M2 M3];
% 
% ddq = InvM*(u-CandM);
% ------------------------------------------------------------------------

%end

% -------------------------------------------------------------------------
function Q = invdyn(dh, rb, qc, qcdot, qcddot, grav, tau, w0, wdot0, a0)
% Inverse dynamic with recursive Newton-Euler following De luca slides
import casadi.*

if nargin < 6
    grav = [0;0;0];
end

z0 = [0; 0; 1];


    q = qc;
 qdot = qcdot;
qddot = qcddot;

w     = SX.sym('w'   , 3,3);
wdot  = SX.sym('wdot', 3,3);
vdot  = SX.sym('vdot', 3,3);
n     = SX.sym('n', 3,3);
F     = SX.sym('F', 3,1);
N     = SX.sym('N', 3,1);

f     = SX.sym('f', 3,3);

Q     = SX.sym('Q', 3,1);

N_DOFS = length(q);
R0 = eye(3);
    % ---------------------------------------------------------------------
    % Forward recursion 
    for i = 1 : N_DOFS
        T = calc_transformation(i-1, i, dh, q);
        R = T(1:3, 1:3);
        %R0 = R0*R;
        p = [dh.a(i); dh.d(i)*sin(dh.alpha(i)); dh.d(i)*cos(dh.alpha(i))];
        z0 = R*z0;
        if i > 1
            w(:,i)    = R'*(w(:,i-1) + z0.*qdot(i));
            wdot(:,i) = R'*(wdot(:,i-1) +  z0.*qddot(i) + ...
                        cross(w(:,i-1), z0.*qdot(i)));
            vdot(i,:) = R'*vdot(:,i-1) + cross(wdot(:,i), p) + ...
                        cross(w(:,i-1), cross(w(:,i),p));
        else
            w(:,i)    = R'*(w0 + z0.*qdot(i));
            wdot(:,i) = R'*(wdot0 + z0.*qddot(i));
            vdot(:,i) = -R'*(a0-grav) + cross(wdot(:,i), p) + ...
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
        z0 = [0;0;1];
        if dh.type(i) == 't'
            Q(i) = f(:,i)'*R'*z0;
        elseif dh.type(i) == 'r'
            Q(i) = n(:,i)'*R'*z0;
        end
    end
    
    f0 = f(:,1);
    n0 = n(:,1);
    
    % Allocation matrix + draft of backward recursion
%     G = [  -0.2449    0.2673   -0.2832    0.2496    0.6017   -0.5546
%             0.4993   -0.4947   -0.4691    0.4971    0.0212   -0.0038
%             0.7157    0.6958    0.6572    0.6714    0.7051    0.7304
%            -0.0031    0.1762    0.1769   -0.0046   -0.1827   -0.1746
%            -0.2100   -0.0996    0.0945    0.2112    0.1111   -0.1077
%             0.1976   -0.1857    0.1864   -0.1945    0.1894   -0.1923    ];
%         
%     fp = m_base*(vdot0-grav)+f0;
%     tp = n0;

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
    if dh.type(i) == 'r'
        dh.theta(i) = qc(i);
    elseif dh.type(i) == 'p'
        dh.d(i) = qc(i);
    end
    
    ct = cos(dh.theta(i) + dh.offset(i));
    st = sin(dh.theta(i) + dh.offset(i));
    ca = cos(dh.alpha(i));
    sa = sin(dh.alpha(i));
    
    T = T * [ ct    -st*ca   st*sa     dh.a(i)*ct ; ...
              st    ct*ca    -ct*sa    dh.a(i)*st ; ...
              0     sa       ca        dh.d(i)    ; ...
              0     0        0         1          ];
end


% %% Initialize
% g = 9.81; % gravity vector
% 
% O3 = zeros(3); % squared zero matrix of order 3
% I3 = eye(3); % identity matrix of order 3
% e3 = [0; 0; 1]; % z unit vector
% 
% % Generalised mass matrix
% M = [m*I3 O3;
%       O3  J ];
% 
% % Allocation matrix
% G = [  -0.2449    0.2673   -0.2832    0.2496    0.6017   -0.5546
%         0.4993   -0.4947   -0.4691    0.4971    0.0212   -0.0038
%         0.7157    0.6958    0.6572    0.6714    0.7051    0.7304
%        -0.0031    0.1762    0.1769   -0.0046   -0.1827   -0.1746
%        -0.2100   -0.0996    0.0945    0.2112    0.1111   -0.1077
%         0.1976   -0.1857    0.1864   -0.1945    0.1894   -0.1923    ];
%  
% % Rotation matrix from Euler Angles
% r11 = cos(theta) * cos(psi);
% r21 = cos(theta) * sin(psi);
% r31 = -sin(theta);
% r12 = sin(phi) * sin(theta) * cos(psi) - cos(phi) * sin(psi);
% r22 = sin(phi) * sin(theta) * sin(psi) + cos(phi) * cos(psi);
% r32 = sin(phi) * cos(theta);
% r13 = cos(phi) * sin(theta) * cos(psi) + sin(phi) * sin(psi);
% r23 = cos(phi) * sin(theta) * sin(psi) - sin(phi) * cos(psi);
% r33 = cos(phi) * cos(theta);
% R = [[r11 r21 r31].', [r12 r22 r32].', [r13 r23 r33].'];  
% 
% % Mapping matrix from eta_dot to omega
% t11 = 1; t12 = 0;         t13 = -sin(theta);
% t21 = 0; t22 = cos(phi);  t23 = sin(phi)*cos(theta);
% t31 = 0; t32 = -sin(phi); t33 = cos(phi)*cos(theta);
% T = [t11 t12 t13; t21 t22 t23; t31 t32 t33];
% 
% %% Dynamic equation of a floating-base rigid-body system
% 
% dyn_eq = M \ ([-m*g*e3; -cross(omega, J*omega)] + [R O3; O3 I3]*G*gamma);
% 
% %% Return output
% % p_dot = p_dot; % unnecessary 
% p_ddot = dyn_eq(1:3);
% eta_dot = T \ omega;
% omega_dot = dyn_eq(4:6);

end
