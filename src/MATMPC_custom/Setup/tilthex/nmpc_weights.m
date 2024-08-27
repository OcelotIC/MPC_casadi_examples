%% Save already-existenting variables
workspaceVars = who; % get worspace variables' names



%% Default weights
weights.q_pos   = [ 10      10      10      ];
weights.q_att   = [ 10      10      10      ];
weights.q_vel   = [ 1       1       1       ];
weights.q_avel  = [ 1       1       1       ];
weights.q_acc   = [ 1e-4    1e-4    1e-4    ];
weights.q_aacc  = [ 1e-4    1e-4    1e-4    ];

%% Collect weights
weights.q = [ weights.q_pos weights.q_att ...
              weights.q_vel weights.q_avel ...
              weights.q_acc weights.q_aacc     ];

%% Export
indata.weights = weights;



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
