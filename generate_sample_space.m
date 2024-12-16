%% This file implements an algorithm for sampling uniformly-in-space the solution to a hierarchical system of nonlinear inequalities
    % The algorithm was originally described in "Mechanisms for dysregulation of excitatory-inhibitory balance underlying allodynia in dorsal horn neural subcircuits"
    % by Ginsberg et al 2024.
% This specific implementation samples uniformly-in-space the system of
    % inequalities defining the allowable parameter space for the "simple"
    % circuit, described in Ginsberg et al 2024.
% To run this particular implementation, the user will need to compute maximum and minimum values for 
    % each coupling strength (g_abi, g_ie, and g_abe) underlying the circuit, or to use the provided values 
    % (see "g_simple_overall_maxs.mat" and "g_simple_overall_mins.mat").

%% Load parameters
master_parameter_file

%%I parameters
v_i_min = beta_i_simple - 12*alpha_i_simple; %the minimum allowed voltage (mV) of the I population for properly functioning simple subcircuit instantiations
    v_i_rest = -60; %the resting voltage (mV) of the inhibitory population
    params_simple(14) = v_i_rest;
    v_i_thr = beta_i_simple - alpha_i_simple; %the "threshold" voltage (mV) of the J population, above which the I pop. fires at a significant level
    f_i_v_i_thr = activation_method(v_i_thr, max_i_simple, alpha_i_simple, beta_i_simple); %firing rate (Hz) of the I pop. at the firing threshold voltage
    v_i_max = beta_i_simple + 12*alpha_i_simple; %the maximum allowed voltage (mV) of the I population for properly functioning simple subcircuit instantiations
    f_i_v_i_max = activation_method(v_i_max, max_i_simple, alpha_i_simple, beta_i_simple); %the firing rate (Hz) of the I population at v_i_max (is nearly 80 Hz)

%%E parameters
v_e_min = beta_e_simple - 12*alpha_e_simple;
    resting_excit_simple = -60;
    v_e_rest = beta_e_simple - 5.4430*alpha_e_simple; %-60
    v_e_thr = beta_e_simple - alpha_e_simple; %-20.9
    v_e_max = beta_e_simple + 12*alpha_e_simple; %77.8

%%f_ab parameters
f_ab_min = 10; %the minimum "typical" average firing rate (Hz) along ABeta fibers in response to non-painful stimuli
f_ab_max = 20; %the maximum "typical" average firing rate (Hz) along ABeta fibers in response to non-painful stimuli
f_ab_params = [f_ab_min; f_ab_max]; %assemble f_ab parameters into one vector

%Assemble parameters into the vectors "e_info" containing excitatory population parameters 
%"i_info" containing inhibitory population parameters and
%"f_ab_info" containing f_ab parameters
    e_info(1) = v_e_min;
    e_info(2) = v_e_rest;
    e_info(3) = v_e_thr;
    e_info(4) = v_e_max;
    i_info(1) = v_i_min;
    i_info(2) = v_i_rest;
    i_info(3) = v_i_thr;
    i_info(4) = v_i_max;
    i_info(5) = alpha_i_simple;
    i_info(6) = beta_i_simple;
    i_info(7) = max_i_simple;
    f_ab_info(1) = f_ab_min;
    f_ab_info(2) = f_ab_max;

g_abi_bounds = [(v_i_thr - v_i_rest)/f_ab_min; (v_i_max - v_i_rest)/f_ab_max];

%% Load the Overall Min and Max Values of each Coupling Strength--these need to be pre-computed.
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


%% Find Points in the Simple Circuit
%Do so by dividing the normalized parameter space into hyperrectangles
%Choose number of samples per unit volume
samples_per_unit_volume = 5000;

%First choose lengths of hyper-rectangles
number_g_abi_intervals = 25;
    g_abi_rectangle_width = 1 / number_g_abi_intervals;
    g_abi_vector = 0 : g_abi_rectangle_width : 1; %the partition of g_abi values

number_g_ie_intervals_per_g_abi = 25;

number_g_abe_intervals = 10;
max_rectangle_volume = 1/number_g_abi_intervals/number_g_ie_intervals_per_g_abi/number_g_abe_intervals;


% Actually do the sampling:
tic %time the algorithm
sample_number = 1;
sample_points_simple = nan*ones(samples_per_unit_volume, 3);
while sample_number <= samples_per_unit_volume
    g_abi_index = ceil(rand*number_g_abi_intervals); %randomly select an index (and therefore a subinterval) in the partition for normalized g_abi values
    g_abi_simple_normalized = g_abi_vector(g_abi_index); %compute the corresponding normalized g_abi value (take the left endpoint of the subinterval in the partition)
        g_abi_simple = g_abi_simple_normalized*(g_abi_overall_max - g_abi_overall_min) + g_abi_overall_min; %compute the corresponding un-normalized g_abi
    
    % given g_abi, compute maximimum and minimum possible g_ie values and compute parameters for a partition of [g_ie_min, g_ie_max]:
    g_ie_min = g_ie_min_fun(g_abi_simple, params_simple, v_e_thr, v_e_rest, f_ab_min); %compute the minimum possible value of g_ie given g_abi
        g_ie_min_normalized = (g_ie_min - g_ie_overall_min)/(g_ie_overall_max - g_ie_overall_min); %normalize that min possible value of g_ie
    g_ie_max = g_ie_max_fun(g_abi_simple, params_simple, e_info, i_info, f_ab_info); %compute the maximimum possible value of g_ie given g_abi
        g_ie_max_normalized = (g_ie_max - g_ie_overall_min)/(g_ie_overall_max - g_ie_overall_min); %normalize that max possible value of g_ie
    g_ie_length = g_ie_max_normalized - g_ie_min_normalized; %compute the length of the interval of possible g_ie values given g_abi
    g_ie_rectangle_width = g_ie_length/number_g_ie_intervals_per_g_abi; %compute the width of subintervals in the partition of possible normalized g_ie values given g_abi
    g_ie_vector = g_ie_min_normalized : g_ie_rectangle_width : g_ie_max_normalized + g_ie_rectangle_width; %the partition of g_ie values
    % randomly select g_ie as follows:
    g_ie_index = ceil(rand*(number_g_ie_intervals_per_g_abi +1) ); %randomly select an index (and therefore a subinterval) in the partition for normalized g_ie
    g_ie_simple_normalized = g_ie_vector(g_ie_index); %compute the corresponding normalized g_ie value (take the left endpoint of the subinterval)
    g_ie_simple = g_ie_simple_normalized*(g_ie_overall_max - g_ie_overall_min) + g_ie_overall_min; %compute the corresponding unnormalized g_ie
        
        %given g_abi and g_ie, compute maximum and minimum possible g_abe values and compute parameters for a partition of [g_ie_min, g_ie_max]:
        g_abe_min = g_abe_min_fun(g_abi_simple, g_ie_simple, e_info, i_info, f_ab_info); %compute the minimum possible value of g_abe given g_abi and g_ie
        g_abe_min_normalized = (g_abe_min - g_abe_overall_min)/(g_abe_overall_max - g_abe_overall_min); %normalize that min possible value of g_abe
        g_abe_max = g_abe_max_fun(g_abi_simple, g_ie_simple, params_simple, v_e_max, v_e_rest, f_ab_max); %compute the maximum possible value of g_abe given g_abi and g_ie
        g_abe_max_normalized = (g_abe_max - g_abe_overall_min)/(g_abe_overall_max - g_abe_overall_min); %normalize that max possible value of g_abe
        g_abe_length = g_abe_max_normalized - g_abe_min_normalized; %compute the length of the interval of possible g_abe values given g_abi and g_ie
        if g_abe_length > 0 %proceed if the interval actually exists given g_abi and g_ie
            g_abe_rectangle_width = g_abe_length/number_g_abe_intervals; %compute the width of subintervals in the partition of possible normalized g_abe values given g_abi and g_ie
            g_abe_vector = g_abe_min_normalized : g_abe_rectangle_width : g_abe_max_normalized + g_abe_rectangle_width; %the partition of g_abe values
            % randomly select g_abe as follows (in the same manner g_abi and g_ie were randomly selected):
            g_abe_index = ceil(rand*(number_g_abe_intervals + 1)); %randomly select an index (and therefore a subinterval) in the partition for g_ie
            g_abe_simple_normalized = g_abe_vector(g_abe_index); %compute the corresponding normalized g_ie value (take the left endpoint of the subinterval)
            g_abe_simple = g_abe_simple_normalized*(g_abe_overall_max - g_abe_overall_min) + g_abe_overall_min; %compute the corresponding unnormalized g_ie
            
            %compute the volume of the hyperrectangle selected
            volume_rectangle = g_abi_rectangle_width*g_ie_rectangle_width*g_abe_rectangle_width;
            number_points_in_rectangle = volume_rectangle/1*samples_per_unit_volume;
            if rand < volume_rectangle/max_rectangle_volume %sample uniformly-at-random from the selected hyper-rectangle with probability proportional to its volume
                g_abi_value = rand*g_abi_rectangle_width + g_abi_simple_normalized; %uniformly-at-random select g_abi from its range
                g_ie_value = rand*g_ie_rectangle_width + g_ie_simple_normalized; %uniformly-at-random select g_ie from its range
                g_abe_value = rand*g_abe_rectangle_width + g_abe_simple_normalized; %uniformly-at-random select g_abe from its range
                if g_ie_min_normalized <= g_ie_value && g_ie_value <= g_ie_max_normalized && g_abe_min_normalized < g_abe_value && g_abe_value < g_abe_max_normalized
                    %if the randomly selected point truly satisfies the inequalities defining the APS, the algorithm proceeds:
                    sample_points_simple(sample_number, :) = [g_abi_value, g_ie_value, g_abe_value]; %store the sampled point
                    sample_number = sample_number + 1; %increment the sample number by 1
                    if mod(sample_number, 100) == 0 %print the sample_number to show how far along the simulation is.
                        sample_number %the simulation is done when sample_number = 5000 is printed
                    end
                end
            end
        end
end
toc

% save the results of the sample as a matrix  called "sample_points_simple"
save('sample_points_simple.mat', 'sample_points_simple')

%Un-normalize points in the sample for future use and save the result as a matrix called "simple_points"
simple_points = sample_points_simple;
simple_points(:, 1) = sample_points_simple(:, 1)*(g_abi_overall_max - g_abi_overall_min) + g_abi_overall_min;
simple_points(:, 2) = sample_points_simple(:, 2)*(g_ie_overall_max - g_ie_overall_min) + g_ie_overall_min;
simple_points(:, 3) = sample_points_simple(:, 3)*(g_abe_overall_max - g_abe_overall_min) + g_abe_overall_min;
save('simple_points.mat', 'simple_points')





















%% Functions
function out = min_slope_f_i(g_abi_simple, params) 
    %have to be from params_simple
    %g_abe_simple = params(1);          % Coupling strength between E population and signals from ABeta fibers
    %g_ce_simple = params(2);           % Coupling strength between E population and signals from C fibers
    %g_ie_simple = params(3);           % Coupling strength between E population and signals sent from I population 
    %g_abi_simple = params(4);          % Coupling strength between I population and signals from ABeta fibers
    %g_nmdac_simple = params(5);        % Coupling strength between NMDA signals arising from C fibers and E pop.
    %alpha_e_simple = params(6);        % Flatness of firing-rate response of E population to changes in E pop.'s voltage
    %beta_e_simple = params(7);         % Average voltage of E population required to make E population fire on average at half its max firing-rate
    alpha_i_simple = params(8);        % Flatness of firing-rate response of I population to changes in I pop.'s voltage
    beta_i_simple = params(9);         % Average voltage of I population required to make E population fire on average at half its max firing-rate
    %max_e_simple = params(10);          % Max firing-rate response of E population to changes in E pop.'s voltage
    %max_i_simple = params(11);         % Max firing-rate response of I population to changes in E pop.'s voltage
    %tau_e_simple = params(12);         % Time constant for E population. Determines the speed at which corresponding diff eqs evolves and responds to stimuli
    %tau_i_simple = params(13);         % Time constant for I population. Determines the speed at which corresponding diff eqs evolves and responds to stimuli
    v_i_rest = params(14); % Steady state voltage of inhibition population in absence of ABeta stimuli
    %resting_exc_simple = params(15);   % Steady state voltage of excitatory population in absence of ABeta stimuli
    %alpha_m_simple = params(16);           %Flatness of response in magnesium deblockage to voltage increases
    %beta_m_simple = params(17);            %average voltage of E population needed to remove half of magnesium blockage
    %tau_m_simple = params(18);             %Time constant for magnesium blockage. Determines the speed at which magnesium blockages is removed/recovers

    f_ab_min = 10;
    f_ab_max = 20;
    x_0 = g_abi_simple*f_ab_min;
    
    lambertw_inputs = -exp(2*(v_i_rest - beta_i_simple)/alpha_i_simple + 1);
    x_1 = alpha_i_simple/2*(1 - lambertw(lambertw_inputs));
        x_1(find(lambertw_inputs<-1/exp(1))) = nan;
        x_1 = real(x_1);
    x_2 = alpha_i_simple/2*(1 - lambertw(-1, -exp(2*(v_i_rest - beta_i_simple)/alpha_i_simple + 1)));
        x_2(find(lambertw_inputs<-1/exp(1))) = nan;
        x_2 = real(x_2);
    x_3 = g_abi_simple*f_ab_max;

    slope_x_0 = activation_method_inhibitory(x_0 + v_i_rest, alpha_i_simple, beta_i_simple)./x_0;
    slope_x_3 = activation_method_inhibitory(x_3 + v_i_rest, alpha_i_simple, beta_i_simple)./x_3;
    min_matrix = min(slope_x_0, slope_x_3);
    size_g_abi_simple = max(size(g_abi_simple));
    if ~isnan(x_1)
        slope_x_1 = activation_method_inhibitory(x_1 + v_i_rest, alpha_i_simple, beta_i_simple)./x_1;
        for g_abi_index_1 = 1 : size_g_abi_simple
            for g_abi_index_2 = 1 : size_g_abi_simple
                x_0_val = x_0(g_abi_index_1, g_abi_index_2);
                x_3_val = x_3(g_abi_index_1, g_abi_index_2);
                if x_0_val < x_1 && x_1 < x_3_val
                    min_matrix(g_abi_index_1, g_abi_index_2) = min(min_matrix(g_abi_index_1, g_abi_index_2), slope_x_1);
                end
            end
        end
    end
    if ~isreal(x_2)
        slope_x_2 = activation_method_inhibitory(x_2 + v_i_rest, alpha_i_simple, beta_i_simple)./x_2;
        for g_abi_index_1 = 1 : size_g_abi_simple
            for g_abi_index_2 = 1 : size_g_abi_simple
                x_0_val = x_0(g_abi_index_1, g_abi_index_2);
                x_3_val = x_3(g_abi_index_1, g_abi_index_2);
                if x_0_val < x_2 && x_2 < x_3_val
                    min_matrix(g_abi_index_1, g_abi_index_2) = min(min_matrix(g_abi_index_1, g_abi_index_2), slope_x_2);
                end
            end
        end
    end
    out = min_matrix;
end


function out = g_ie_min_fun(g_abi_simple, params_simple, v_e_lthr, v_e_rest, f_ab_min)
    out = (v_e_lthr - v_e_rest)/f_ab_min/min_slope_f_i(g_abi_simple, params_simple)/g_abi_simple;
end

function g_ie_max = g_ie_max_fun(g_abi_simple, params_simple, e_info, i_info, f_ab_info)
    v_e_min = e_info(1);
    v_e_rest = e_info(2);
    v_e_lthr = e_info(3);
    v_e_max = e_info(4);
    %v_i_min = i_info(1);
    v_i_rest = i_info(2);
    %v_i_lthr = i_info(3);
    %v_i_max = i_info(4);
    alpha_i_simple = i_info(5);
    beta_i_simple = i_info(6);
    max_i_simple = i_info(7);
    f_ab_min = f_ab_info(1);
    f_ab_max = f_ab_info(2);
    
    v_i_max = g_abi_simple*f_ab_max + v_i_rest;
    f_i_max = activation_method(v_i_max, max_i_simple, alpha_i_simple, beta_i_simple);
    
    g_ie_min = g_ie_min_fun(g_abi_simple, params_simple, v_e_lthr, v_e_rest, f_ab_min);
    g_ie_max_zeroing_function = @(g_ie_simple)(g_abe_max_fun(g_abi_simple, g_ie_simple, params_simple, v_e_max, v_e_rest, f_ab_max) - g_abe_min_fun(g_abi_simple, g_ie_simple, e_info, i_info, f_ab_info));
    g_ie_max = fzero(g_ie_max_zeroing_function, [g_ie_min + 0.0000001, g_ie_min + 100]);
    g_ie_max = max(g_ie_max, (v_e_rest - v_e_min)/f_i_max);
end




function out = g_abe_max_fun(g_abi_simple, g_ie_simple, params_simple, v_e_max, v_e_rest, f_ab_max)
    out = min(g_abi_simple*g_ie_simple*min_slope_f_i(g_abi_simple, params_simple), (v_e_max - v_e_rest)/f_ab_max);
end

function out = g_abe_min_fun(g_abi_simple, g_ie_simple, e_info, i_info, f_ab_info)
    v_e_min = e_info(1);
    v_e_rest = e_info(2);
    v_e_lthr = e_info(3);
    %v_e_max = e_info(4);
    %v_i_min = i_info(1);
    v_i_rest = i_info(2);
    %v_i_lthr = i_info(3);
    %v_i_max = i_info(4);
    alpha_i_simple = i_info(5);
    beta_i_simple = i_info(6);
    max_i_simple = i_info(7);
    f_ab_min = f_ab_info(1);
    f_ab_max = f_ab_info(2);
    
    [~, out] = extrema_of_bad_function(g_abi_simple, g_ie_simple, alpha_i_simple, beta_i_simple, max_i_simple, v_i_rest, -(v_e_rest - v_e_min), f_ab_min, f_ab_max);
    
    g_abe_min = (v_e_lthr - v_e_rest)/f_ab_min;
    out = max(out, g_abe_min);
end




