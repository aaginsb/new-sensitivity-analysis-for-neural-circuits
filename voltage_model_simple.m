function dxdt = voltage_model_simple(t, x, t_vec_ab, avg_fr_ab, t_vec_c, avg_fr_c, params)
% Simple circuit
% Inputs: 
    % t = times of interest (t),
    % x = voltages of any populations plus mag if applicable (see below)-- (x)
    % t_vec_ab, avg_fr_ab = times, average number of spikes in corresponding time intervals incoming from ABeta fibers
    % t_vec_c, avg_fr_c = times, average number of spikes in corresponding time intervals incoming from C fibers. avg_fr_c = 0 in simulations from the paper
    % params--constants used in the model--see below
% Outputs: time derivatives of v_e, v_i at times in t

%% Create shorthand for variables
v_e = x(1); % Average voltage of excitatory population
v_i = x(2); % Average voltage of inhibitory population
mag = x(3); % Average proportion of full magnesium blockage on E cells due to C inputs

%% Vector of Parameters 
g_abe = params(1);          % Coupling strength between E population and signals from ABeta fibers
g_ce = params(2);           % Coupling strength between E population and signals from C fibers
g_ie = params(3);           % Coupling strength between E population and signals sent from I population 
g_abi = params(4);          % Coupling strength between I population and signals from ABeta fibers
g_nmdac = params(5);        % Coupling strength between NMDA signals arising from C fibers and E pop.
%alpha_e = params(6);        % Flatness of firing-rate response of E population to changes in E pop.'s voltage
%beta_e = params(7);         % Average voltage of E population required to make E population fire on average at half its max firing-rate
alpha_i = params(8);        % Flatness of firing-rate response of I population to changes in I pop.'s voltage
beta_i = params(9);         % Average voltage of I population required to make E population fire on average at half its max firing-rate
%max_e = params(10);          % Max firing-rate response of E population to changes in E pop.'s voltage
max_i = params(11);         % Max firing-rate response of I population to changes in E pop.'s voltage
tau_e = 1.2*params(13); %params(12);         % Time constant for E population. Determines the speed at which corresponding diff eqs evolves and responds to stimuli
tau_i = params(13);         % Time constant for I population. Determines the speed at which corresponding diff eqs evolves and responds to stimuli
resting_inhib = params(14); % Steady state voltage of inhibition population in absence of ABeta stimuli
resting_exc = params(15);   % Steady state voltage of excitatory population in absence of ABeta stimuli
alpha_m = params(16);           %Flatness of response in magnesium deblockage to voltage increases
beta_m = params(17);            %Average voltage of E population needed to remove half of magnesium blockage
tau_m = params(18);             %Time constant for magnesium blockage. Determines the speed at which magnesium blockages is removed/recovers


%% time dependent input functions - do not change
f_ab = interp1(t_vec_ab, avg_fr_ab, t); % Interpolate the data set (ft, f) at times t
f_c = interp1(t_vec_c, avg_fr_c, t); % Interpolate the data set (gt, g) at times t

%% response curves 
f_i = activation_method(v_i, max_i, alpha_i, beta_i); %avg firing rates of tonic inhibitory population.
m_inf = activation_method(v_e, 1.0, alpha_m, beta_m); %Minf

%% differential equations for VOLTAGE
d_ve = (g_abe*f_ab + (g_ce + g_nmdac*mag)*f_c - g_ie*f_i - v_e + resting_exc)/ tau_e;  % E
d_vi = (g_abi*f_ab - v_i + resting_inhib) / tau_i; % I
d_mag = (m_inf - mag) / tau_m;  % magnesium--the default in simulations is to set f_c = 0, making mag irrelevant in those simulations

dxdt = [d_ve; d_vi; d_mag];

end
        
