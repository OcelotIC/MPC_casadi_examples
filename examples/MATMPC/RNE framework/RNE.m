
clc
clear all
close all

addpath('/home/ichelikh/rvctools/robot')
addpath('/home/ichelikh/rvctools/common') 
addpath('/home/ichelikh/rvctools/smtb')

figure;
hold on;
h1 = plot3(0,0,0, 'b');   % The link
h2 = plot3(0,0,0, 'or');  % The joint/end-effector
h3 = plot3(0,0,0, 'xk');  % The center of mass
xlabel('x')
ylabel('y')
zlabel('z')
view([1 1 1]);
axis equal;
xlim([-1 1]);
ylim([-1 1]);
zlim([-1 1]);
legend('Link', 'Joint', 'Center of Mass')

% -------------------------------------------------------------------------
% Define the DH parameters
N_DOFS = 3;
dh.theta = [0 0 0];
dh.alpha = [pi/2 0 0];
dh.offset = [pi 0 0];
dh.d = [0 0 0];
dh.a = [0 0.3 0.3];
dh.type = ['r' 'r' 'r'];

% -------------------------------------------------------------------------
% Rigid body paramaters: inertia, mass, and cener of mass
I1 = diag([0.0005  0.000011 0.0004]);
I2 = diag([0.00035 0.0018   0.0016]);
I3 = diag([0.423   1.414    1.185])*1e-4;
rb.I =  repmat(eye(3), 1, 1, N_DOFS);

rb.I(:,:,1) = I1;
rb.I(:,:,2) = I2;
rb.I(:,:,3) = I3;
rb.m = [0.059 0.256 0.093];
% In standard DH, COM is mesured respect to its end-effector (using local
% frame). When COM is [0 0 0]', it means the COM is located exactly at the
% end-effector. Therefore, COM usually has negative values, which means the
% COM is behind the end-effector
rb.r = [0 0 0; -0.15 0 0; -0.15 0 0]';

% -------------------------------------------------------------------------
% Arbitrary trajectory as the inputs: joint position, velocity, and 
% acceleration
ts = 0.001;
time_span = 0:ts:1;
qc = [ pi/4*sin(4*pi*1*time_span)' 0*time_span' pi/4*sin(4*pi*1*time_span)']; %pi/6*sin(2*pi*1*time_span)'
qcdot = gradient(qc', ts)';
qcddot = gradient(qcdot', ts)';

% -------------------------------------------------------------------------
% Kinematic simulation, optional, for visualization purpose!   
% End-effector position, form the base which is located at [0 0 0]'
EE = zeros(3, N_DOFS+1);
COM = zeros(3, N_DOFS);    
for k = 1 : length(time_span)  
    for i = 1 : 1 : N_DOFS
        T = calc_transformation(0, i, dh, qc(k,:));
        EE(:,i+1) = T(1:3,4);
        COM(:,i) = EE(:,i+1) + T(1:3,1:3) * rb.r(:,i);
    end
    x1 = EE(1, 4);
    y1 = EE(2, 4);
    z1 = EE(3, 4);
    [x1 y1 z1]
    % Draw the robot
    set(h1, 'XData', EE(1, :), 'YData', EE(2, :),'ZData', EE(3, :));
    set(h2, 'XData', EE(1, :), 'YData', EE(2, :),'ZData', EE(3, :));
    set(h3, 'XData', COM(1, :), 'YData', COM(2, :),'ZData', COM(3, :));
    drawnow;
end

% -------------------------------------------------------------------------
% Here we go!
Q = invdyn(dh, rb, qc, qcdot, qcddot, [0; 0; -9.8]);

% ------------------------------------------------------------------------
% Plotting
figure;
hold on;
plot(time_span, Q(1,:), 'b');
plot(time_span, Q(2,:), 'r');
plot(time_span, Q(3,:), 'g');

% -------------------------------------------------------------------------
% For validation purpose, using RVC Toolbox for comparison
% Disable these lines below if you don't want to try Peter Corkee's RVC 
% toolbox
L(1) = Revolute('d', pi, 'a', 0, 'alpha', pi/2, ...
    'm', rb.m(1), 'I', rb.I(:,:,1), 'r',  rb.r(:,1));
L(2) = Revolute('d', 0, 'a', 0.3, 'alpha', 0, ...
    'm', rb.m(2), 'I', rb.I(:,:,2), 'r', rb.r(:,2));
L(3) = Revolute('d', 0, 'a', 0.3, 'alpha', 0, ...
    'm', rb.m(3), 'I', rb.I(:,:,3), 'r', rb.r(:,3));

threelink = SerialLink(L, 'name', 'two link');
threelink.gravity = [0; 0; -9.8];
torque_by_rvc_toolbox = threelink.rne(qc, qcdot, qcddot);
plot(time_span, torque_by_rvc_toolbox(:,1), '-.b', 'LineWidth', 2);
plot(time_span, torque_by_rvc_toolbox(:,2), '-.r', 'LineWidth', 2);
plot(time_span, torque_by_rvc_toolbox(:,3), '-.g', 'LineWidth', 2);

legend('Torque 1', 'Torque 2', 'Torque 3', ...
     'Torque 1 by RVC', 'Torque 2 by RVC', 'Torque 3 by RVC');


% -------------------------------------------------------------------------
function Q = invdyn(dh, rb, qc, qcdot, qcddot, grav)
% Inverse dynamic with recursive Newton-Euler

if nargin < 6
    grav = [0;0;0];
end

z0 = [0; 0; 1];

for k = 1 : length(qc)  
    q = qc(k, :);
    qdot = qcdot(k, :);
    qddot = qcddot(k, :);

    N_DOFS = length(q);
    
    % ---------------------------------------------------------------------
    % Forward recursion
    for i = 1 : N_DOFS
        T = calc_transformation(i-1, i, dh, q);
        R = T(1:3, 1:3);
        p = [dh.a(i); dh.d(i)*sin(dh.alpha(i)); dh.d(i)*cos(dh.alpha(i))];
        
        if i > 1
            w(:, i) =  R'*(w(:, i-1) + z0.*qdot(i));
            wdot(:, i) = R'*(wdot(:, i-1) +  z0.*qddot(i) + ...
                cross(w(:, i-1), z0.*qdot(i)));
            vdot(:,i) = R'*vdot(:,i-1) + cross(wdot(:,i), p) + ...
                cross(w(:,i), cross(w(:,i),p));
        else
            w(:, i) =  R'*(z0.*qdot(i));
            wdot(:, i) = R'*(z0.*qddot(i));
            vdot(:,i) = R'*grav + cross(wdot(:,i), p) + ...
                cross(w(:,i), cross(w(:,i),p));
        end
    end
    
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
        
        if dh.type(i) == 't'
            Q(i,k) = f(:,i)'*R'*z0;
        elseif dh.type(i) == 'r'
            Q(i,k) = n(:,i)'*R'*z0;
        end
    end
end
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

end