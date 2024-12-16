%% Master Parameter File: weights between populations vary from microcircuit to microcircuit

%% Parameters for the simulation:
% number of realizations:
    number_random_seeds= 1;
% final time (in seconds):
    t_fin = 1.2; %was 0.9
% vector of times at which we output solutions to ODEs being solved
    tspan = 0 : 0.001 : t_fin - 0.001;  % in seconds
% suffix for the name of the output file(s)
    name = '_testNewModel';


%% Parameters for the incoming spike trains on the fibers: (these params are the same for all microcircuits)
% in Hz - baseline FR w.o stimulus
    baseline_freq = 1; 
% in Hz - responses to stimuli
    stim_freq_ab = 20; % was 20;  % response to stimulus on AB fibers
    stim_freq_c = 20; % response to stimulus on C fibers
%Number fibers in bundles
    number_fibers_ab = 300; % was 300; % number of AB fibers in bundle--had been 300
    number_fibers_c = 820; % number of AD and C fibers in bundle
%Times over which stimuli occur
    %ABeta stimuli:
        duration_of_stimulus_ab = 500; % was 500% Time (in ms) during which stimulus is present 
        initial_baseline_duration_ab = 200; % was 200% Time (in ms) passed before which stimulus arrives
    %C stimuli:
        duration_of_stimulus_c = 500; %Was 500 Time (in ms) during which stimulus is present 
        initial_baseline_duration_c = 350;  %Was 500 Time (in ms) during which stimulus is present 

% put all this information into a vector to be passed into a function file:
    stim_info_ab = [baseline_freq, stim_freq_ab, number_fibers_ab, initial_baseline_duration_ab, duration_of_stimulus_ab];
    stim_info_c = [baseline_freq, stim_freq_c, number_fibers_c, initial_baseline_duration_c, duration_of_stimulus_c];


%% Parameters for the basic circuit
% Vector of Parameters 
g_abe_simple = 1;              %1 Coupling strength between E population and signals from ABeta fibers
g_ce_simple = 3.0;               %2 Coupling strength between E population and signals from C fibers
g_ie_simple = 1;              %3 Coupling strength between E population and signals sent from I population 
g_abi_simple = 1;              %4 Coupling strength between I population and signals from ABeta fibers
g_nmdac_simple = 3;        %5 Coupling strength between NMDA signals arising from C fibers and E pop.
alpha_e_simple = 7.9;         %6 Flatness of firing-rate response of E population to changes in E pop.'s voltage
beta_e_simple = -17;           %7 Average voltage of E population required to make E population fire on average at half its max firing-rate
alpha_i_simple = 9.3;          %8 Flatness of firing-rate response of I population to changes in I pop.'s voltage
beta_i_simple = -30;           %9 Average voltage of I population required to make E population fire on average at half its max firing-rate
max_e_simple = 50.0;           %10 Max firing-rate response of E population to changes in E pop.'s voltage
max_i_simple = 80.0;             %11 Max firing-rate response of I population to changes in E pop.'s voltage--has identical effect on excitatory voltage as rescaling e
tau_e_simple = 0.01;             %12 Time constant for E population. Determines the speed at which corresponding diff eqs evolves and responds to stimuli
tau_i_simple = .02;             %13 Time constant for I population. Determines the speed at which corresponding diff eqs evolves and responds to stimuli
resting_inhib_simple = -60;    %14 Steady state voltage of inhibition population in absence of ABeta stimuli--choosing this so the basically firing I population is still firing, but only a little in the absence of ABeta input
resting_exc_simple = -60; %15 Steady state voltage of excitatory population in absence of ABeta stimuli, was -70
alpha_m_simple = 15.2; %16
beta_m_simple = -40.0; %17 %was originally -20, set to -40 because voltage was never coming close to -20
tau_m_simple = 3; %18
    
params_names_simple = {'g_{abe}--simple' 'g_{ce}--simple' 'g_{ie}--simple' 'g_{abi}--simple' 'g_{nmdac}--simple' 'alpha_{e}--simple' 'beta_{e}--simple' 'alpha_{i}--simple' 'beta_{i}--simple' 'max_{e}--simple' 'max_{i}--simple' 'tau_{e}--simple' 'tau_{i}--simple' 'resting inhib--simple' 'resting exc--simple' 'alpha_m--simple' 'beta_m--simple' 'tau_m--simple'};
params_simple = [g_abe_simple g_ce_simple g_ie_simple g_abi_simple g_nmdac_simple alpha_e_simple beta_e_simple alpha_i_simple beta_i_simple max_e_simple max_i_simple tau_e_simple tau_i_simple resting_inhib_simple resting_exc_simple alpha_m_simple beta_m_simple tau_m_simple];

% initial voltage for each population
ve_0_simple = resting_exc_simple; %should be fixed at -70
vi_0_simple = resting_inhib_simple; %should be fixed at resting_inhib_simple. Doing so would make v_I start at resting_inhib_simple and return to resting_inhib_simple in absence of ABeta input
mag_0_simple = 0;
x_0_simple = [ve_0_simple; vi_0_simple; mag_0_simple];



