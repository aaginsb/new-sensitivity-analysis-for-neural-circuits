%% Data Analysis
% This file will plot the uniform-in-space random sample of points, (computed using "generate_sample_space.m"), 
    % from the allowable parameter space for the simple circuit (see "Mechanisms for dysregulation of
    % excitatory-inhibitory balance underlying allodynia in dorsal horn neural subcircuits" by Ginsberg et al 2024).
% This file will also plot time courses of the voltages and firing-rates for the excitatory population 
    % and the inhibitory population in the particular instantiations of the simple subcircuit.


%% Parameters

clear

master_parameter_file

% f_ab parameters
f_ab_min = 10;
f_ab_max = 20;



%% Load the sample points
load('sample_points_simple.mat')
sample_points = sample_points_simple;

load('simple_points.mat')

%Load the Overall Min and Max Values of each Coupling Strength
load('g_simple_overall_mins.mat')
gs_overall_mins = g_simple_overall_mins;
g_abi_overall_min = g_simple_overall_mins(1);
g_ie_overall_min = g_simple_overall_mins(2);
g_abe_overall_min = g_simple_overall_mins(3);

load('g_simple_overall_maxs.mat')
gs_overall_maxs = g_simple_overall_maxs;
g_abi_overall_max = g_simple_overall_maxs(1);
g_ie_overall_max = g_simple_overall_maxs(2);
g_abe_overall_max = g_simple_overall_maxs(3);

sample_count = length(simple_points(:, 1));

% Set the default tick label to use latex formatting
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')







%% Analyze Data from Raw Sample Points
%% Plot the Sample Space
figure
scatter3(sample_points(:, 1), sample_points(:, 2), sample_points(:, 3))
xlabel('$\hat{g}_{ABI}$', 'interpreter', 'latex')
ylabel('$\hat{g}_{IE}$', 'interpreter', 'latex')
zlabel('$\hat{g}_{ABE}$', 'interpreter', 'latex')
xlim([0, 1])
ylim([0, 1])
zlim([0, 1])
set(gca, 'fontsize', 20)


%% Do an example simulation of the model
% Choose which point from the sample to examine more closely
sample_number = 1;

% Extract the coupling strengths corresponding to the current sample, and
% reorder them for use in "voltage_model_simple"
params_test = simple_points(sample_number, :);
params_simple([1, 3, 4]) = params_test([3, 2, 1]);

%% Generate a noisy ABeta stimulus to drive this particular instantiation of the simple circuit

% Set parameters for generating the noisy Abeta stimulus
stim_freq_ab = rand*(f_ab_max - f_ab_min) + f_ab_min; %choose a random Abeta firing rate in the normal range of ABeta firing rates
stim_freq_c = 20; %the default value of stim_freq_c--we'll set the firing rate along the c fibers to 0 shortly
stim_freq_c_end = 20; %we'll set the firing rate along the c fibers to 0 shortly
stim_info_ab(2) = stim_freq_ab;
stim_info_ab(5) = 500; %make simulation last 50 ms

% actually get the firing rate
[t_vec_ab, avg_fr_ab] = make_smooth_avg_fr(stim_info_ab, t_fin); %average firing rate for ABeta fibers--mechanoreception. Same for all microcircuits
[t_vec_c, avg_fr_c] = make_smooth_avg_fr(stim_info_c, t_fin); %average firing rate for C fibers--nociception. Same for all MCs
avg_fr_c = 0*avg_fr_c;


% Plot the average firing rate of ABeta fires
figure
plot(t_vec_ab, avg_fr_ab, 'color', 'k', 'linewidth', 2)
xlabel('Time (s)', 'interpreter', 'latex')
ylabel('$f_{A\beta}$ (Hz)', 'interpreter', 'latex')
hold off
set(gca, 'fontsize', 20)



%% Simulate the simple system, driving it with the average firing-rates of the ABeta fires simulated above.
[t_simple,x_simple] = ode45(@(t_simple,x_simple) voltage_model_simple(t_simple, x_simple, t_vec_ab, avg_fr_ab, t_vec_c, avg_fr_c, params_simple), tspan, x_0_simple);






%% Plot the voltages and firing rates of the inhibitory and excitatory populations

%plot the voltage of the inhibitory population
figure
plot(t_simple, x_simple(:, 2), 'color', 'r', 'linewidth', 2, 'linewidth', 2)
xlabel('Time (s)', 'interpreter', 'latex')
ylabel('$V_{I}$ (mV)', 'interpreter', 'latex')
hold off
set(gca, 'fontsize', 20)

% Plot the firing rate of the inhibitory population
firing_rate = activation_method(x_simple(:, 2), max_i_simple, alpha_i_simple, beta_i_simple);     
figure
plot(t_simple, firing_rate, 'color', 'r', 'linewidth', 2)
xlabel('Time (s)', 'interpreter', 'latex')
ylabel('$f_{I}$ (Hz)', 'interpreter', 'latex')
hold off
set(gca, 'fontsize', 20)

% Plot the voltage of the excitatory population
figure
plot(t_simple, x_simple(:, 1), 'color', 'b', 'linewidth', 2)
xlabel('Time (s)', 'interpreter', 'latex')
ylabel('$V_{E}$ (mV)', 'interpreter', 'latex')
hold off
set(gca, 'fontsize', 20)

% Plot the firing rate of the excitatory population
firing_rate = activation_method(x_simple(:, 1), max_e_simple, alpha_e_simple, beta_e_simple);     
figure
plot(t_simple, firing_rate, 'color', 'b', 'linewidth', 2)
xlabel('Time (s)', 'interpreter', 'latex')
ylabel('$f_{E}$ (Hz)', 'interpreter', 'latex')
hold off
set(gca, 'fontsize', 20)
