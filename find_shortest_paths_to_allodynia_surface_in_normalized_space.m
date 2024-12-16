%% Find the Shortest Path to Allodynia Surface in Normalized Parameter Space


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load parameters, sampled points from the allowable parameter space (APS) for the
% simple circuit, and the maximum and minimum of coupling strengths used to
% normalize the coupling strengths

clear
master_parameter_file %load all parameters for the simple circuit. 
    % Note that master_parameter_file may contain some parameters unnecessary for this file.

% I parameters
    v_i_rest = -60; % resting voltage of the inhibitory population
    v_i_lthr = beta_i_simple - alpha_i_simple; % voltage of the inhibitory population above which
        % a significant firing rate is assumed to occur. 
    v_i_max = beta_i_simple + 12*alpha_i_simple; % maximum allowable voltage for inhibitory
        % populations in properly function circuits

% E parameters
    v_e_rest = -60; % resting voltage of the excitatory population
    v_e_lthr = beta_e_simple - alpha_e_simple; % voltage of the excitatory population above which
        % a significant firing rate is assumed to occur. 

% f_ab parameters
f_ab_min = 10;
f_ab_max = 20;



%Load the sample points
load('sample_points_simple.mat')
sample_points = sample_points_simple;

load('simple_points.mat')



% Load the Overall Min and Max Values of each Coupling Strength
load('g_simple_overall_mins.mat')
g_abi_overall_min = g_simple_overall_mins(1);
g_ie_overall_min = g_simple_overall_mins(2);
g_abe_overall_min = g_simple_overall_mins(3);

load('g_simple_overall_maxs.mat')
g_abi_overall_max = g_simple_overall_maxs(1);
g_ie_overall_max = g_simple_overall_maxs(2);
g_abe_overall_max = g_simple_overall_maxs(3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Actually find the shortest paths from the sampled allowable parameter space (APS)  points to the allodynia surface
number_sampled_points = length(sample_points(:, 1));

% Initialize an array to store the nearest points on the allodynia surface
    %corresponding to each sampled point from the allowable parameter space:
nearest_point_on_allodynia_surface_data_normalized = nan*ones(number_sampled_points, 3);

% Initialize an array to store the displacements from the allowable parameter space points to the 
    % corresponding nearest points on the allodynia surface:
displacement_to_allodynia_surface_data_normalized = nan*ones(number_sampled_points, 3);

%Initialize an array to store the corresponding distances.
distance_to_allodynia_surface_data_normalized = nan*ones(number_sampled_points, 1);

tic
for simple_point_index = 1 : number_sampled_points
    % access the current APS point:
    point_normalized = sample_points_simple(simple_point_index, :).';
        % g_abi = point_normalized(1);
        % g_ie = point_normalized(2);
        % g_abe = point_normalized(3);

    %find the nearest point on the allodynia surface
    nearest_point_on_allodynia_surface = shortest_normalized_distance_to_allodynia_surface_multi_start(point_normalized, v_i_rest, alpha_i_simple, beta_i_simple, v_e_rest, v_e_lthr, alpha_e_simple, beta_e_simple, g_abi_overall_min, g_ie_overall_min, g_abe_overall_min, g_abi_overall_max, g_ie_overall_max, g_abe_overall_max, f_ab_min, f_ab_max);
        % g_abi = nearest_point_on_allodynia_surface(1);
        % g_ie = nearest_point_on_allodynia_surface(2);
        % g_abe = nearest_point_on_allodynia_surface(3);
        % f_ab = nearest_point_on_allodynia_surface(4);

    %store the nearest point on the allodynia surface:
    nearest_point_on_allodynia_surface_data_normalized(simple_point_index, :) = nearest_point_on_allodynia_surface(1:3);
    %store the displacement to the nearest point on the allodynia surface:
    displacement_to_allodynia_surface_data_normalized(simple_point_index, :) = nearest_point_on_allodynia_surface(1:3) - point_normalized;
    %store the distance to the nearest point on the allodynia surface:
    distance_to_allodynia_surface_data_normalized(simple_point_index) = norm(nearest_point_on_allodynia_surface(1:3) - point_normalized);
end
toc
% Save the nearest points on the allodynia surface, the corresponding 
    % displacements from the allowable parameter space points to the nearest
    % points on the allodynia surface, and the corresponding distances
save('nearest_point_on_allodynia_surface_data_normalized.mat', 'nearest_point_on_allodynia_surface_data_normalized')
save('displacement_to_allodynia_surface_data_normalized.mat', 'displacement_to_allodynia_surface_data_normalized')
save('distance_to_allodynia_surface_data_normalized.mat', 'distance_to_allodynia_surface_data_normalized')

% Use a parallel plot to show 500 of the shortest paths (i.e., displacements), 
    % from APS points to their nearest points on the allodynia surface
figure
for sample_number = 1 : 500
    plot([1, 2, 3], displacement_to_allodynia_surface_data_normalized(sample_number, :), 'color', 'k');
    hold on
end
hold off
xlim([.9, 3.1])
ylim([-0.5, 0.5])
ylabel('$\Delta$(Normalized Coupling Strengths)', 'interpreter', 'latex')
xticks(linspace(1, 3, 3))
xticklabels({'$\hat{g}_{A\beta I}$', '$\hat{g}_{IE}$', '$\hat{g}_{A\beta E}$'})
set(gca, 'fontsize', 20)































%% Functions
function [true_X_new, objective_true_X_new] = shortest_normalized_distance_to_allodynia_surface_multi_start(point_normalized, v_i_rest, alpha_i_simple, beta_i_simple, v_e_rest, v_e_lthr, alpha_e_simple, beta_e_simple, g_abi_overall_min, g_ie_overall_min, g_abe_overall_min, g_abi_overall_max, g_ie_overall_max, g_abe_overall_max, f_ab_min, f_ab_max) 
    % This function implements the algorithm used to find the shortest path between
        % a point in the allowable parameter space (APS) and the allodynia surface for the simple circuit model,
        % as described in "Mechanisms for dysregulation of excitatory-inhibitory balance 
        % underlying allodynia in dorsal horn neural subcircuits" by Ginsberg et al 2024
    % Inputs to this function:
        % The main input is some point in the allowable parameter space 
            % called "point_normalized."
        % Other inputs such as v_i_rest, alpha_i_simple, etc., are constants.
    % outputs of this function:
        % true_X_new = the [g_abi, g_ie, g_abe, f_ab] vector where (g_abi, g_ie, g_abe) 
            % is the nearest point on the allodynia surface (in normalized parameter space),
            % and f_ab is the normalized ABeta firing-rate for which the simple
            % circuit instantiation with coupling strengths g_abi, g_ie, and g_abe 
            % produces allodynia.
        % (optional): objective_true_X_new = the distance between the APS point and the
            % nearest point on the allodynia surface.
    
    % The objective function below is the squared euclidean distance
        % between the normalized aps point "point_normalized" and a point in
        % space--the first three components of "g_s_and_f_ab".
    objective = @(g_s_and_f_ab) sum((g_s_and_f_ab(1:3) - point_normalized).^2);
        % g_abi = g_s_and_f_ab(1);
        % g_ie = g_s_and_f_ab(2);
        % g_abe = g_s_and_f_ab(3);
        % f_ab = g_s_and_f_ab(4);

    % We will minimize the objective function subject to the constraints that:
        %1. g_abi, and g_ie lie between 2 and their value in 
            % normalized parameter space corresponding to an un-normalized
            % value of 0
        %2. f_ab lies between 0 and 1 (treating f_ab as if it has been normalized)
        %3. g_abe lies between 2 the value of the point_normalized's g_abe
        %4. The E population's voltage is as small as possible so that 
            % E fires at a significant level given f_ab: 
            % v_e_lthr = v_e at steady state = g_abe*f_ab - g_ie*f_i + v_e_rest.
            % This condition is re-written for convenience as 
            % 0= g_abe*f_ab - g_ie*f_i - (v_e_lthr  - v_e_rest).

    % To conduct this minimization, we will Matlab's "fmincon" function,
        % which uses a gradient descent method to find a local minimum.
    % To find the global minimum, we:
        %1. choose a random initial point in the hypercube from which fmincon will start its
            % search, and run fmincon
        %2. Repeat for an additional 15 randomly selected initial points and
            % record which point returned by fmincon outputs the smallest objective
            % function.

    % Define parameters for use in "fmincon"
    X0 = zeros(4, 1); %initial value from which to start the search
    A = []; %there are no linear inequality constraints
    B = []; %there are no linear inequality constraints
    Aeq = []; %there are no linear equality constraints
    Beq = []; %there are no linear equality constraints
    UB = [2*ones(3, 1); 1]; %make the upper bound of the g's 2 
        %and the upper bound of f_ab be 1 (in normalized parameter space)
    LB = [-g_abi_overall_min/(g_abi_overall_max - g_abi_overall_min); -g_ie_overall_min/(g_ie_overall_max - g_ie_overall_min); ...
          -g_abe_overall_min/(g_abe_overall_max - g_abe_overall_min); 0]; %make the minimum of all the g's be 0 and f_ab to be f_ab_min
    LB(3) = point_normalized(3); %set lower bound on g_abe to the value of the point_normalized's g_abe.
        %NONLCON_comp_1 = []; %no nonlinear inequality constraints
    
    % run the local optimization algorithm fmincon with the previously mentioned constraints
    options = optimset('TolCon',1e-10); %set a tolerance criteria to be smaller than the default of 10^(-6)
    X = fmincon(@(g_s_norm) objective(g_s_norm),X0,A,B,Aeq,Beq,LB,UB, @(g_s_norm)...
        NONLCON_norm(g_s_norm, alpha_i_simple, beta_i_simple, v_i_rest, v_e_rest, v_e_lthr, g_abi_overall_max, g_ie_overall_max, g_abe_overall_max, g_abi_overall_min, g_ie_overall_min, g_abe_overall_min, f_ab_min, f_ab_max),...
        options);
    
    % Compute the constraint error. If the constraint error is too
            % large, throw out the solution by replacing the solution X 
            % to the local optimization problem with an arbitrary, large value.
    [~, constraint_error] = NONLCON_norm(X, alpha_i_simple, beta_i_simple, v_i_rest, v_e_rest, v_e_lthr, g_abi_overall_max, g_ie_overall_max, g_abe_overall_max, g_abi_overall_min, g_ie_overall_min, g_abe_overall_min, f_ab_min, f_ab_max);
    if abs(constraint_error)>10^(-4)
        X = 1000*ones(3, 1);
    end
    true_X_new = X; %record the solution to the local optimization problem
    objective_true_X_new = objective(X); %record the value of the objective function given
        %the local mininum X.
        
    number_initial_points_used = 15;
    for corner_number = 1 : number_initial_points_used
        X0 = rand(4,1).*(UB - LB) + LB; %randomly choose a point to start at
            %normalized values of g_abi, g_ie, g_abe, and f_ab
            %will be selected uniformly-at-random between their respective
            %upper and lower bounds.
        
        % run the local optimization algorithm fmincon with the previously mentioned constraints
        options = optimset('TolCon',1e-10);
        X = fmincon(@(g_s_norm) objective(g_s_norm),X0,A,B,Aeq,Beq,LB,UB,...
            @(g_s_norm) NONLCON_norm(g_s_norm, alpha_i_simple, beta_i_simple, v_i_rest, v_e_rest, v_e_lthr, g_abi_overall_max, g_ie_overall_max, g_abe_overall_max, g_abi_overall_min, g_ie_overall_min, g_abe_overall_min, f_ab_min, f_ab_max),...
            options);
        
        % Compute the constraint error. If the constraint error is too
            % large, throw out the solution by replacing the solution X to 
            % the local optimization problem with an arbitrary, large value.
        [~, constraint_error] = NONLCON_norm(X, alpha_i_simple, beta_i_simple, v_i_rest, v_e_rest, v_e_lthr, g_abi_overall_max, g_ie_overall_max, g_abe_overall_max, g_abi_overall_min, g_ie_overall_min, g_abe_overall_min, f_ab_min, f_ab_max);
        if abs(constraint_error)>10^(-4)
            X = 1000*ones(4, 1);
        end

        % Compute the value of the objective function. If the objective
        % function for the solution X found by the local optimization algorithm
        % objective_X is smaller than the previously smallest objective function
        % objective_true_X_new for local solution true_X_new,
        % replace true_X_new with X and objective_true_x_new with objective_X.
        objective_X = objective(X);
        if objective_X < objective_true_X_new && abs(constraint_error) < 10^(-4)
            objective_true_X_new = objective_X;
            true_X_new = X; 
        end
    end
end

function [c, c_eq] = NONLCON_norm(g_s_norm, alpha_i_simple, beta_i_simple, v_i_rest, v_e_rest, v_e_lthr, g_abi_overall_max, g_ie_overall_max, g_abe_overall_max, g_abi_overall_min, g_ie_overall_min, g_abe_overall_min, f_ab_min, f_ab_max)
    %the output of this function are the nonlinear constraints used by 
    %the built-in Matlab function "fmincon" when called in 
    %the function defined above, "shortest_normalized_distance_to_allodynia_surface_multi_start"

    c = [];
    g_s = g_s_norm;
        % g_abi = g_s(1); (where g_abi has been normalized)
        % g_ie = g_s(2); (where g_ie has been normalized)
        % g_abe = g_s(3); (where g_abe has been normalized)
        % f_ab = g_s(4); (where f_ab has been normalized to be in [0, 1]).

    %un-normalize the g_s
    g_s(1) = g_s_norm(1)*(g_abi_overall_max - g_abi_overall_min) + g_abi_overall_min;
    g_s(2) = g_s_norm(2)*(g_ie_overall_max - g_ie_overall_min) + g_ie_overall_min;
    g_s(3) = g_s_norm(3)*(g_abe_overall_max - g_abe_overall_min) + g_abe_overall_min;
    g_s(4) = g_s_norm(4)*(f_ab_max - f_ab_min) + f_ab_min; %normalize f_ab as well
    
    % The function "shortest_normalized_distance_to_allodynia_surface_multi_start"
        % minimizes an objective function subject to the constraint that c_eq = 0, where
        % c_eq = g_abe*f_ab - g_ie*f_i - (v_e_lthr  - v_e_rest).
    % This equation for c_eq is simply the condition that the E population's
        % voltage is as small as possible so that E fires at a significant level
        % given f_ab: v_e_lthr = v_e at steady state = g_abe*f_ab - g_ie*f_i + v_e_rest.
    c_eq = g_s(3)*g_s(4) - g_s(2)*activation_method_inhibitory(g_s(1)*g_s(4) + v_i_rest, alpha_i_simple, beta_i_simple)...
        - (v_e_lthr - v_e_rest);   
end

