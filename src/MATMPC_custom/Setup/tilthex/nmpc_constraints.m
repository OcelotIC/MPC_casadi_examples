%% Save already-existenting variables
workspaceVars = who; % get worspace variables' names



%% Load from Settings
nc = settings.nc;        % Number of constraints
ncN = settings.ncN;      % Number of constraints at terminal stage
nu = settings.nu;        % Number of inputs

%% Load from data
% Propeller coefficients
c_f = indata.model.c_f;

%% Specify state bounds
% Bounds on forces
f_min = c_f * 16^2;
f_max = c_f * 102^2;

% +20% tollerance in f_min and f_max bounds (for actuators health)
f_min = f_min; % + 0.2 * f_min;
f_max = f_max; % - 0.2 * f_max;

% Bounds on position and attitude
x_min = -5;
x_max = 5;
y_min = -2;
y_max = 2;
z_min = 0.3;
z_max = 5;
rp_min = -pi/2; % roll and pitch
rp_max = pi/2;

%% Collect state bounds
lb_x = [x_min; y_min; z_min; rp_min; rp_min; f_min * ones(nu, 1)];
ub_x = [x_max; y_max; z_max; rp_max; rp_max; f_max * ones(nu, 1)];

%% Input bounds (here only params; computation in Simulink)
% Each element of omega_bar is related to one value of 
% dot_omega_bar_max and one of dot_omega_bar_min. 
omega_bar = [30, 40, 50, 60, 70, 80, 90]';  % Hz

% dot_omega_bar_min=[-60,-80,-100,-70,-80,-80,-70]';
dot_omega_bar_min = [-127, -121, -114, -118, -128, -111, -95]';

% dot_omega_bar_max=[100,100,100,80,90,90,90]';
dot_omega_bar_max = [209, 208, 244, 208, 149, 156, 135]';

% These coefficients (coeff_angular_acc and _dec) are used for 
% linear interpolation of values of omega which are in-between 
% of two values inside omega_bar vector.
coeff_angular_dec = zeros(6,1);
coeff_angular_acc = zeros(6,1);

for i = 1:1:6  
    coeff_angular_dec(i) = (dot_omega_bar_min(i+1) - dot_omega_bar_min(i)) / (omega_bar(i+1) - omega_bar(i));
    coeff_angular_acc(i) = (dot_omega_bar_max(i+1) - dot_omega_bar_max(i)) / (omega_bar(i+1) - omega_bar(i));
end


dot_alpha_min = deg2rad(-100); % -10 [deg/s]
dot_alpha_max = deg2rad(100);  %  10 [deg/s]

%% Export
indata.constraints.lbx = lb_x;
indata.constraints.ubx = ub_x;
indata.constraints.index.z = 3; % Used to dynamically relax z lower bound at the beginning of real fight
indata.constraints.omega_bar = omega_bar;
indata.constraints.dot_omega_bar_min = dot_omega_bar_min;
indata.constraints.dot_omega_bar_max = dot_omega_bar_max;
indata.constraints.coeff_angular_dec = coeff_angular_dec;
indata.constraints.coeff_angular_acc = coeff_angular_acc;
indata.constraints.f_min = f_min;
indata.constraints.f_max = f_max;

indata.constraints.dot_alpha_min = dot_alpha_min;
indata.constraints.dot_alpha_max = dot_alpha_max;

%% Specify variables of this script to keep
vars2keep = {'indata'};

%% Workspace cleanup
for i= 1:length(vars2keep)
  if ~any(strcmp(workspaceVars, vars2keep{i}))
    workspaceVars = union(workspaceVars, vars2keep{i});
  end
end
variables2remove = setdiff(who, workspaceVars);
clear(variables2remove{:}); % clean all unecessary variables created here
clear variables2remove
