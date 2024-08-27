%% Save already-existenting variables
workspaceVars = who; % get worspace variables' names


%% Define the list of waypoints with corresponding average speed

% --- Positions of waypoints (w.r.t. World Frame)
%wps.pos = [0;0;3];
wps.pos = [ 0  1 0;% 1 0 0];
            1  1 1;%0 1 1]; 
            3  3 3];%1 2 3] ];

%WPs = [WPs [    0;  -0.8; 2.1]];
wps.q   = [pi 5*pi/6  pi;%,   pi,  pi,  pi];
           0  pi/3    0 ;%, 0.5, 0.2,  0] ;
           0  pi/4    0];%, 0.3, 0.1,  0] ];
       
%wps.q= [pi;0;0];
wps.num = size(wps.pos, 2); % number of columns of array WPs

% --- Speeds between each waypoints
wps.vel = 1*ones(1, wps.num);

% Time constants
wps.t_start = 0.1;    % Beginning of the first motion
wps.t_rest = 1;     % Pause between two motions

%% Export
indata.wps = wps;

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