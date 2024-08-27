close all;
clear all;
clc;

load('urdf.mat'); % load condensed urdf

%% README
% CONDENSED URDF vs CONTROLLER MODEL ======================================
%
% The difference between the two is the different frame definitions. It may
% be the cause that in the urdf you have the same number of frames as in
% the robot model, but they might be differently placed and oriented.
%
% PROCEDURE ===============================================================
%
% The purpose of this script is to convert the inertial parameters
% (inertia tensor and com location*) of the condensed urdf into those
% needed by the robot model.
% The main idea is to exploit the fact that the body frame is coincident
% and expressed in the same way for the condensed urdf and the controller
% model. In fact, it is assumed that both are oriented in the same way and
% placed in the com of the base** (UAV).
% So the procedure is to take all the necessary inertial parameters from
% the condensed urdf and expressed them in its body frame.
% Since the two body frame are coincident, these parameters are coeherent
% also in the robot model.
% At this point, exploiting the knowledge of how the frames are expressed
% in the robot model, it is possible to rotate back the inertial
% parameters of each link in its own frame of the robot model.
%
% * the mass in uneffected.
% ** if necessary, this last assumption can be removed. If so, it is necessary
%    to express properly the inertia tensor.
% =========================================================================

%% DH parameters of robot model
% DOC =====================================================================
%
% These parameters will position and orientation of each frame in the robot 
% model.
% =========================================================================
a1 = 0; % length of link 1
a2 = 0.31; % length of link 2
a3 = 0.21; % length of link 3
% Values of the DH parameters when the robot is in its home (rest)
% configuration
dh.a = [a1 a2 a3];
dh.d = [0 0 0];
dh.alpha = [pi/2 0 0];
dh.theta = [pi 0 0];

A0 = [ eye(3)      zeros(3,1); % B_A_0 : from Link Frame 0 to Base Frame
       zeros(1,3)      1      ];
A = HTML2L(dh.a, dh.d, dh.alpha, dh.theta);
A1 = A(1:4,:); % 0_A_1 : from Link Frame 1 to Link Frame 0
A2 = A(5:8,:); % 1_A_2
A3 = A(9:12,:); % 2_A_3

%% Compute Homogenous Transformation Matrices (robot model)
B_A_l1 = A0*A1;
B_A_l2 = B_A_l1*A2;
B_A_l3 = B_A_l2*A3;

l1_A_B = invHomoMat(B_A_l1);
l2_A_B = invHomoMat(B_A_l2);
l3_A_B = invHomoMat(B_A_l3);

%% Express inertial parameters of the condensed URDF wrt Base Frame
% NOTATION ================================================================
%
% uav: uav frame
% l1: link1 frame
% l2: link2 frame
% l3: link3 frame
% a_T_b: homogenous transformation matrix from frame b to frame a
% =========================================================================
uav_T_l1 = homoMat(att2RotMat('XYZ', urdf.l1.fr_rpy), urdf.l1.fr_xyz);
l1_T_l2 = homoMat(att2RotMat('XYZ', urdf.l2.fr_rpy), urdf.l2.fr_xyz);
l2_T_l3 = homoMat(att2RotMat('XYZ', urdf.l3.fr_rpy), urdf.l3.fr_xyz);

uav_T_l2 = uav_T_l1 * l1_T_l2;
uav_T_l3 = uav_T_l2 * l2_T_l3;

% Express the inertia tensor and the com location of every link wrt Body
% frame of the condensed urdf.
% Those coincides with the inertial parameters of the robot model expressed
% in body frame.
uav = urdf.uav;
l1.B_I = getRotHomoMat(uav_T_l1) * urdf.l1.I * getRotHomoMat(uav_T_l1)';
l1.B_com_xyz = uav_T_l1 * [urdf.l1.com_xyz; 1];

l2.B_I = getRotHomoMat(uav_T_l2) * urdf.l2.I * getRotHomoMat(uav_T_l2)';
l2.B_com_xyz = uav_T_l2 * [urdf.l2.com_xyz; 1];

l3.B_I = getRotHomoMat(uav_T_l3) * urdf.l3.I * getRotHomoMat(uav_T_l3)';
l3.B_com_xyz = uav_T_l3 * [urdf.l3.com_xyz; 1];

%% Move (back) inertial parameters of the robot model from body frame to each link frame
% NOTATION ================================================================
%
% -- Masses do not need any action.
% -- fr_xyz,R in the robot model are the position and orientation of each 
%    link wrt the previous one. Thus these information are obtained from
%    the DH parameters.
% -- com_xyz and I are the only inertial parameters that have to be rotated
%    (back) from Body frame to each corresponding link frame.
% =========================================================================
model.uav = urdf.uav;
model.uav.B_A_0 = A0;

model.l1.m = urdf.l1.m; % masses are unaffected!
model.l1.I = getRotHomoMat(l1_A_B) * l1.B_I * getRotHomoMat(l1_A_B)';
model.l1.com_xyz = l1_A_B * l1.B_com_xyz;
model.l1.fr_xyz = getPosHomoMat(A1); % wrt previous link (thus frame Link 0).
model.l1.fr_R = getRotHomoMat(A1);

model.l2.m = urdf.l2.m;
model.l2.I = getRotHomoMat(l2_A_B) * l2.B_I * getRotHomoMat(l2_A_B)';
model.l2.com_xyz = l2_A_B * l2.B_com_xyz;
model.l2.fr_xyz = getPosHomoMat(A2);
model.l2.fr_R = getRotHomoMat(A2);

model.l3.m = urdf.l3.m;
model.l3.I = getRotHomoMat(l3_A_B) * l3.B_I * getRotHomoMat(l3_A_B)';
model.l3.com_xyz = l3_A_B * l3.B_com_xyz;
model.l3.fr_xyz = getPosHomoMat(A3);
model.l3.fr_R = getRotHomoMat(A3);

%% Save final model
save('model.mat', 'model');
clear all;