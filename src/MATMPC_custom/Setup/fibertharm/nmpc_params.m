%% Save already-existenting variables
workspaceVars = who; % get worspace variables' names



%% Export
settings.N = 15;
settings.N2 = 1;  % 5   % No. of partial condensing blocks --- Used if set opt.condesing = 'partial_condensing' 
settings.r = 1;   % 10  % No. of input move blocks --- Used if opt.qpsolver = 'qpoases_mb' (mb: moving blocks)

%% Options
settings.opt.integrator = 'ERK4-CASADI'; % 'ERK4','IRK3, 'ERK4-CASADI'
settings.opt.condensing = 'default_full'; % 'default_full','no','blasfeo_full(require blasfeo installed)','partial_condensing'
settings.opt.qpsolver = 'qpoases'; 
settings.opt.hotstart = 'no'; % 'yes','no' (only for qpoases)
settings.opt.shifting = 'no'; % 'yes','no'
settings.opt.lin_obj = 'no'; % 'yes','no' % if the inner objective function is linear and the outer objective is sum of quadratic
settings.opt.ref_type = 1; % 0-time invariant, 1-time varying(no preview), 2-time varying (preview)
settings.opt.nonuniform_grid = 0; % supports only ERK4 and IRK3
settings.opt.enable.outputVars = 1; % Enables computation of output variables inside NMPC block

%% Specify variables of this script to keep
vars2keep = {'settings'};



%% Workspace cleanup
for i= 1:length(vars2keep)
  if ~any(strcmp(workspaceVars, vars2keep{i}))
    workspaceVars = union(workspaceVars, vars2keep{i});
  end
end
variables2remove = setdiff(who, workspaceVars);
clear(variables2remove{:}); % clean all unecessary variables created here
clear variables2remove