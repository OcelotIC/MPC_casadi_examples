
%% Setup

clc;
clear;
close all;
rng(0);

%% Crawling torso problem definition

% Initial state (q_b, q, q_b_dot, q_dot)
x_0 = [0.0; 0.0; 0.0; 0.0];

% Target state 
x_star = [0.0; 0.0; pi; 0.0];

% Time horizon
t_f = 2.0;

% Number of states along trajectory
N = floor(t_f ./ 0.01);

% Maximum magnitude of control
u_max = 0.1*size(tau);

% Initialize dynamics
fprintf("initializing crawl_vispa dynamics...\n")
m_c = 1.0;
m_p = 0.01;
l = 0.25;
dyn = WalkingDynamics(robot);

% Initialize cost
fprintf("initializing quadratic cost function...\n")
% Q_f = [10.0 0 0 0;
%        0 10.0 0 0;
%        0 0 10.0 0;
%        0 0 0 10.0];
q_b_weights = 10*ones(1,6);
q_weights = 10*ones(1,robot.n_q);
q_b_dot_weights = 10*ones(1,6);
q_dot_weights = 10*ones(1,robot.n_q); 
%lambdas_weights = 10*ones(1,6);
Q_f = diag([q_b_weights, q_weights, q_b_dot_weights, q_dot_weights]);
%Q_f = diag([q_b_weights, q_weights, q_b_dot_weights, q_dot_weights, lambdas_weights]);
R = 0.01*ones(1, robot.n_q);
cost = QuadraticCost(Q_f, R);

% Number of DDP iterations
num_iter = 1000;

% DDP learning rate
alpha = 0.1;

% Video framerate
fps = 30;

%% Execution of DDP

fprintf("executing DDP...");

tic;
sol = ddp(x_0, x_star, t_f, N, dyn, cost, u_max, num_iter, alpha);
toc;

%% Begin post-processing of solution

if sol.error == 1
    fprintf("DIVERGENCE ERROR: try decreasing learning rate\n");
    return;
end

% Extract the state from the solution
z = zeros(1, length(sol.x));
theta = zeros(1, length(sol.x));
z_dot = zeros(1, length(sol.x));
theta_dot = zeros(1, length(sol.x));
u = zeros(1, length(sol.x));
for k = 1:N
    z(k) = sol.x{k}(1);
    theta(k) = sol.x{k}(3);
    z_dot(k) = sol.x{k}(2);
    theta_dot(k) = sol.x{k}(4);
    u(k) = sol.u{k};
end
