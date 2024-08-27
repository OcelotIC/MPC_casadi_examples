% ------------------------------------------------------
%  This is an example of initializing simulink simulation
%  ------------------------------------------------------

%%
clear mex; close all; clc;

addpath([pwd,'/nmpc']);
addpath([pwd,'/model_src']);
addpath([pwd,'/mex_core']);
addpath([pwd,'/solver/win64/qpoases']);
addpath([pwd,'/robotics']);
%% Parametri Simulazione
cd data;
if exist('settings','file')==2
    load('settings');
    load('indata');
    %Gdrone = load('alloc_mat_fibertharm_urdf');
    cd ..
else 
    cd ..
    error('No setting data is detected!');
end
%disp(Gdrone)
%missing parameters in new model gen for fibertharm these are added here
%but will be removed asap
Ts = 0.015;

nu_lambda = indata.model.nu_lambda ;
nu_v = indata.model.nu_v;

settings.Ts = Ts;
settings.nu_lambda = nu_lambda;
settings.nu_v = nu_v;

%Ts  = settings.Ts;    % Sampling time


Ts_st = settings.Ts_st;  % Shooting interval
nx = settings.nx;        % No. of differential states
nu = settings.nu;        % No. of controls
nz = settings.nz;        % No. of algebraic states
ny = settings.ny;        % No. of outputs (references)    
nyN= settings.nyN;       % No. of outputs at terminal stage 
np = settings.np;        % No. of parameters (on-line data)
nc = settings.nc;        % No. of constraints
ncN = settings.ncN;      % No. of constraints at terminal stage
nbu = settings.nbu;      % No. of control bounds
nbx =   settings.nbx;      % No. of state bounds
nbu_idx = settings.nbu_idx; % Index of control bounds
nbx_idx = settings.nbx_idx; % Index of state bounds

%% add more to Settings

N  = 15;
N2 = 10; 
r  = 10;

settings.N = N;
settings.N2 = N2;
settings.r = r;

%% options
opt.hessian   ='Gauss_Newton';  % 'Gauss_Newton', 'Generalized_Gauss_Newton'
opt.integrator='ERK4'; % 'ERK4','IRK3, 'IRK3-DAE'
opt.condensing='default_full';  %'default_full','no','blasfeo_full(require blasfeo installed)','partial_condensing'
opt.qpsolver  ='qpoases'; 
opt.hotstart  ='no'; %'yes','no' (only for qpoases)
opt.shifting  ='no'; % 'yes','no'
opt.ref_type  = 0; % 0-time invariant, 1-time varying(no preview), 2-time varying (preview)
opt.nonuniform_grid= 0; % currently not supported 
opt.RTI            = 'yes'; % if use Real-time Iteration

%% available qpsolver
%'qpoases' (for full condensing)
%'qpoases_mb' (for full condensing+moving block, please use ERK4 as the integrator)
%'hpipm_sparse' (run mex_core/compile_hpipm.m first; set opt.condensing='no')
%'hpipm_pcond' (run mex_core/compile_hpipm.m first; set opt.condensing='no')
 
%% Initialization
%must change this either with a outside function like Gianluca or other
%% Load from Settings
nc = settings.nc;        % Number of constraints
ncN = settings.ncN;      % Number of constraints at terminal stage
nu = settings.nu;        % Number of inputs
np = settings.np;

%% Load from data
% Propeller coefficients
c_f = indata.model.c_f;
% Mass and gravity intensity

if strcmp(model, 'tilthex')
    m = indata.m;
else
    m = indata.model.rb.mdrone;
end
g = 9.81; %indata.model.g;
% First waypoint

cd Setup;

if strcmp(model, 'tilthex')
    cd tilthex;

else   
    
    cd fibertharm;
end

%run('nmpc_params'); % useless
run('nmpc_constraints');
run('nmpc_weights');
run('nmpc_wps');
run('nmpc_start');
cd ..
cd ..
wp0 = indata.wps.pos(:,1);

%% 

pos_t0 = wp0;
vel_t0 = zeros(3, 1);
att_t0 = [0; 0; 0];
omega_t0 = zeros(3, 1);

if strcmp(model, 'fibertharm')
    q0  = indata.wps.q(:,1);
    dq0 = zeros(3,1);
end

%compute gravity force and torque at the base with respect of CoM offset

% 
% p_CoM = compute_CoM(pos_t0, att_t0, q0)
% 
% init_force = -m*g*[0;0;1;  p_CoM];
% gamma_t0 = Gdrone.G\init_force;

hover_thrust = m * g / (1 * cos(20*pi/180) * 6);
gamma_t0 = hover_thrust * ones(nu_lambda, 1);

%% Collect
% Initial state vector
if strcmp(model, 'fibertharm')
    x0 = [pos_t0; vel_t0; att_t0; omega_t0; q0; dq0];
else
    x0 = [pos_t0; vel_t0; att_t0; omega_t0; gamma_t0];
end
% Initial input vector
u0 = zeros(nu, 1);
u0(1:6) = hover_thrust;
% Initial parameter vector
para0 = zeros(max(1, np), 1);

%% Export
indata.x0 = x0;
indata.u0 = u0;
indata.para0 = para0;

%W=repmat([10 10 0.1 0.1 0.01]',1,N);
%WN=W(1,1:nyN); %switched indexes need to talk to this with Gianluca

%% constraints 
% upper and lower bounds for states (=nbx)


lb_x = indata.constraints.lbx;
ub_x = indata.constraints.ubx ;

               
 % upper and lower bounds for general constraints (=nc)
 lb_g = [];
 ub_g = [];            
 lb_gN = [];
 ub_gN = [];

%%
 lb = repmat(lb_g,N,1);
 ub = repmat(ub_g,N,1);
 lb = [lb;lb_gN];
 ub = [ub;ub_gN];

%% Define vectors' sizes --- Needed by Simulink
if (N * nc + ncN) > 0
   settings.size_b = [(N * nc + ncN)  1]; % Size of the vector for general constraints
else
    settings.size_b = 1; % scalar
end

settings.size_bu = [nu N]; % Size of the vector for input constraints

if nbx > 0
    % settings.size_bx = [nbx N+1]; % Size of the vector for state constraints
    settings.size_bx = [nbx N];
else
    % settings.size_bx = 1; % scalar
    settings.size_bx = [1 N];
end

if np > 0
    settings.size_od = [np N+1];
else
    settings.size_od = [1 N+1];
end

settings.size_W = [ny N];
settings.size_WN = nyN;
settings.size_ref = [ny N];
settings.size_refN = nyN;
settings.size_x = [nx N];
settings.size_xN = [nx 1];

opt.enable_outputVars = 1;
z0 = zeros(nz,1);
z = repmat(z0,1,N);
if isempty(z)
    z0=0;
    z=0;
end
settings.z0 = z0;