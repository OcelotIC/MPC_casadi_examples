%script for executing everything if ur feelin lazy
clear all;
close all;
clc

% addpath('C:\Users\idris\ms2021-idriss_chelikh\Software\robotics')
%addpath('/home/ichelikh/ms2021-idriss_chelikh/Software/robotics')
model = 'fibertharm'; 

if strcmp(model, 'tilthex')
    simulink_model = 'MATMPC_SIMULINK_TILTHEX';
else 
    simulink_model = 'MATMPC_SIMULINK_FIBERTHARM';
end

Initialization_Simulink;

disp(' ');
fprintf('** Run %s.slx\n  Please wait, this may take a while...\n', ...
            simulink_model);

settings.sim_time = indata.wps.num * 5;
tic
sim(simulink_model, settings.sim_time);
toc

% Check if simulation has been completed
if time >= settings.sim_time
    disp('* Simulation Completed!');
    outdata.sim_complete = true;
else
    time
    disp('* Simulation terminated before reaching the end!');
    outdata.sim_complete = false; 
end

%% Parse data
% disp(' ');
% fprintf('** Parse Simulink output...\n');
% parse_outdata;
% disp('* done');
plot_data
fprintf('\nEnd!\n');