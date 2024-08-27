close all;
clear all;
clc;

%% README
% ROBOT'S LINKS ===========================================================
%
%   root_link <--> arm_base <--> arm_link1 <--> arm_link2 <--> arm_link3
%
% Each link its connected to its predecessor, on the left, through a joint
% which defines the position and orientation of that link's frame wrt its
% predecessor (URDF convention).
% 
% URDF CONVENTION =========================================================
%
% -- JOINT
% The joint frame is coincident with the child-link frame.
% The xyz and rpy defines the position and orientation of the joint
% (child-link) frame wrt the parent one.
%
% -- LINK
% The xyz and rpy fields in the inertial, visual and collision sub-elements
% are defined wrt the link frame (which is defined by a joint which has
% this link as child in its declaration).
%
% -------------------------------------------------------------------------
% Therefore each link will have its own frame (link frame), whose position 
% and orientation wrt  the previous one, is specified by joint, and an 
% inertial frame (link-inertial frame) defined by the fields xyz and rpy in
% the sub-element inertial (similarly for collisions and visuals).
%
% So the frame parameters are taken from the joint declarations, while the
% inertial parameters from the link ones (see URDF).
%
% DOC =====================================================================
%
% <link_name> : in {root_link, arm_base, arm_link1, arm_link2, arm_link3}
%            .fr : frame parameters
%               .xyz = position of the link frame's origin
%               .rpy = orientation of the link frame
%            .in : inertial parameters
%               .m = mass
%               .I = inertia tensor expressed in inertial-link frame
%               .xyz = origin of the link-inertial frame expressed in link
%                      frame
%               .rpy = orientation of the link-inertial frame wrt link
%                      frame (it provides the rotation matrix link_R_linkInertial)
% =========================================================================

%% Parameters from original URDF
% Uav 
root_link.fr.xyz = [0; 0; 0]; % wrt World Frame
root_link.fr.rpy = [0; 0; 0];
root_link.in.m = 2.08;
root_link.in.I = diag([0.0115, 0.0114, 0.0194]);
root_link.in.xyz = [0; 0; 0];
root_link.in.rpy = [0; 0; 0];

% Base of the 3DoF manipulator
arm_base.fr.xyz = [0; 0; -0.1046]; % wrt root_link
arm_base.fr.rpy = [-1.5708; -1.5708; -2.0944];

arm_base.in.xyz = [0; 0; -0.04981];  % pos of arm_base's com wrt arm_base's frame
arm_base.in.rpy = [0; 0; 0]; % relative orientation between inertial frame and link frame
arm_base.in.m = 0.512;
arm_base.in.I = [0.00089561525  0.000000033869084  -0.000000050910484;
                 0.000000033869084  0.00089680894  0.000033914855;
                -0.000000050910484 0.000033914855 0.00081659326];

% Link 1
arm_link1.fr.xyz = [0; 0; 0]; % wrt arm_base
arm_link1.fr.rpy = [0; 0; 0];

arm_link1.in.xyz = [0; 0; 0]; 
arm_link1.in.rpy = [1.5707; 0; 0];
arm_link1.in.m = 0.059;
arm_link1.in.I = [0.000048053428  -0.0000000021805677  0.0000000033360114;
                 -0.0000000021805677  0.000010977023  -0.000000011384151;
                  0.0000000033360114 -0.000000011384151 0.000039181274];

% Link2
arm_link2.fr.xyz = [0; 0; 0]; % wrt arm_link1
arm_link2.fr.rpy = [-1.5708; 3.1416; 0];

arm_link2.in.xyz = [0.155; -0.010; -0.002]; 
arm_link2.in.rpy = [-0.002; 0; 0];
arm_link2.in.m = 0.256;
arm_link2.in.I = [0.00035178136  -0.000036788259  0.0000078669302;
                 -0.000036788259  0.000010977023  -0.0000098704892;
                  0.0000078669302 -0.0000098704892 0.0016336677];

% Link3
arm_link3.fr.xyz = [0.2955; 0; 0]; % wrt arm_link2
arm_link3.fr.rpy = [0; 0; 0];

arm_link3.in.xyz = [0.04; 0; 0]; 
arm_link3.in.rpy = [0; -0.05; 0];
arm_link3.in.m = 0.093;
arm_link3.in.I = [4.2318882e-05  -7.1584145e-08  5.3969447e-06;
                 -7.1584145e-08  1.4139416e-04  2.3843598e-08;
                  5.3969447e-06 2.3843598e-08 1.1853674e-4];


%% Homogenous transformation matrices
% NOTATION ================================================================
%
% rl: root_link
% ab: arm_base
% al1: arm_link1
% al2: arm_link2
% al3: arm_link3
% COM: Center Of Mass
% a_T_b : homogenous transformation matrix from frame b to frame a
% =========================================================================

rl_T_ab = homoMat(att2RotMat('XYZ', arm_base.fr.rpy), arm_base.fr.xyz);

ab_T_al1 = homoMat(att2RotMat('XYZ', arm_link1.fr.rpy), arm_link1.fr.xyz);

al1_T_al2 = homoMat(att2RotMat('XYZ', arm_link2.fr.rpy), arm_link2.fr.xyz);          

al2_T_al3 = homoMat(att2RotMat('XYZ', arm_link3.fr.rpy), arm_link3.fr.xyz);

rl_T_abCOM = rl_T_ab * ...
              homoMat(att2RotMat('XYZ', arm_base.in.rpy), arm_base.in.xyz);

l1_T_l1COM = homoMat(att2RotMat('XYZ', arm_link1.in.rpy), arm_link1.in.xyz);

l2_T_l2COM = homoMat(att2RotMat('XYZ', arm_link2.in.rpy), arm_link2.in.xyz);

l3_T_l3COM = homoMat(att2RotMat('XYZ', arm_link3.in.rpy), arm_link3.in.xyz);

%% Compacting into main rigid bodies, and condensing parameters
% PROCEDURE ===============================================================
%
% Here, it condenses multiple rigid bodies into a single one. Thus, for a
% condesed rigid body its equivalent mass, inertia tensor and com location
% have to be computed.
% -- Equivalent mass = sum of all masses
% -- Equivalent inertia tensor: 1) move all inertia tensors in the same
%                                  frame
%                               2) sum them together
% -- Equivalent com location: weighted-sum of the products between every
%                             rigid body mass and its com location, divided 
%                             by the equivalent mass
% =========================================================================
% UAV = root_link + 3d_arm_base
urdf.uav.fr_xyz = root_link.fr.xyz;
urdf.uav.fr_rpy = root_link.fr.rpy;
urdf.uav.m = root_link.in.m + arm_base.in.m;
urdf.uav.com_xyz = (root_link.in.xyz*root_link.in.m + getPosHomoMat(rl_T_abCOM)*arm_base.in.m)/urdf.uav.m;
urdf.uav.I = root_link.in.I + ...
              getRotHomoMat(rl_T_abCOM) * arm_base.in.I * getRotHomoMat(rl_T_abCOM)';

% Link1 = Link1
urdf.l1.fr_xyz = arm_link1.fr.xyz; % wrt uav
urdf.l1.fr_rpy = arm_link1.fr.rpy;
urdf.l1.m = arm_link1.in.m;
urdf.l1.com_xyz = getPosHomoMat(l1_T_l1COM); % wrt its own frame
urdf.l1.I = getRotHomoMat(l1_T_l1COM) * arm_link1.in.I * getRotHomoMat(l1_T_l1COM)'; % wrt its own frame

% Link2 = Link2
urdf.l2.fr_xyz = arm_link2.fr.xyz; % wrt l1
urdf.l2.fr_rpy = arm_link2.fr.rpy;
urdf.l2.m = arm_link2.in.m;
urdf.l2.com_xyz = getPosHomoMat(l2_T_l2COM); % wrt its own frame
urdf.l2.I = getRotHomoMat(l2_T_l2COM) * arm_link2.in.I * getRotHomoMat(l2_T_l2COM)';

% Link3 = Link3
urdf.l3.fr_xyz = arm_link3.fr.xyz; % wrt l2
urdf.l3.fr_rpy = arm_link3.fr.rpy;
urdf.l3.m = arm_link3.in.m;
urdf.l3.com_xyz = getPosHomoMat(l3_T_l3COM); % wrt its own frame
urdf.l3.I = getRotHomoMat(l3_T_l3COM) * arm_link3.in.I * getRotHomoMat(l3_T_l3COM)';

%% Save condensed URDF
% URDF vs CONDENSED URDF 
% It's called condensed, because it joins together the parameters of some 
% rigid bodies.
%
% Here:
% -- Original URDF: UAV(root_link) + Manipulator Base(arm_base) + link1 + link2 + link3
% -- Condensed URDF: UAV* + link1 + link2 + link3
%
%   where UAV* = root_link + arm_base
% =========================================================================
save('urdf.mat', 'urdf');
clear all;
