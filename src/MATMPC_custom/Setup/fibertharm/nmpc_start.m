%% Save already-existenting variables
workspaceVars = who; % get worspace variables' names



%% Load from Settings
nc = settings.nc;        % Number of constraints
ncN = settings.ncN;      % Number of constraints at terminal stage
nu = settings.nu;        % Number of inputs
np = settings.np;

%% Load from data
% Propeller coefficients
c_f = indata.model.c_f;
% Mass and gravity intensity
m = indata.model.rb.mdrone;
g = 9.81;%indata.model.g;
% First waypoint
wp0 = indata.wps.pos(:,1);

%% 

pos_t0 = wp0;
vel_t0 = zeros(3, 1);
att_t0 = [0; 0; 0];
omega_t0 = zeros(3, 1);

hover_thrust = m * g / (cos(35*pi/180) * cos(10*pi/180) * 6);
gamma_t0 = hover_thrust * ones(nu, 1);

%% Collect
% Initial state vector
x0 = [pos_t0; vel_t0; att_t0; omega_t0; gamma_t0];
% Initial input vector
u0 = zeros(nu, 1);
% Initial parameter vector
para0 = zeros(max(1, np), 1);

%% Export
indata.x0 = x0;
indata.u0 = u0;
indata.para0 = para0;



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
