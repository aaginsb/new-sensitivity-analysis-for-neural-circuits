
function [g_abe_min, g_abe_max] = extrema_of_bad_function(g_abi, g_ie, alpha_i, beta_i, max_i, v_i_rest, C, f_ab_min, f_ab_max)

%% This function finds the extrema of (a*tanh(x-b) + c)/x, where
    % a = 0.5*g_ie*alpha_i/g_abi
    % b = - (v_i,rest - beta_i)/alpha_i
    % c = C*alpha_i/g_abi + a,
% and the particular parameter values could be for example
    % g_abi = 1;
    % g_ie = 1;
    % alpha_i = 9.3;
    % beta_i = -35;
    % v_i_rest = -60;
    % C = 0.01;
    % f_ab_min = 0.01;
    % f_ab_max = 1;


    %compute simpler variables:
    g_abi_divided_by_alpha_i = g_abi/alpha_i;
    a = 0.5*max_i*g_ie*g_abi_divided_by_alpha_i;
    b = -(v_i_rest - beta_i)/alpha_i;
    c = C*g_abi_divided_by_alpha_i + a;
    d = c/a;
    
    x_min = f_ab_min*g_abi_divided_by_alpha_i;
    x_max = f_ab_max*g_abi_divided_by_alpha_i;
    u_min =  x_min - b;
    u_max = x_max - b;
    search_points = [u_min];
    
    %identify extrema of q' (where q is defined as an auxiliary function, below), and check to see if they are in [u_min,
    %u_max]. If so, record them as ``search point''.
    value = -2*b + d - 1;
    lambert_w_input = -2*(1 + d)*exp(2*value);
    if lambert_w_input>=0
        u_0 = real(0.5*value - 0.25*lambertw(lambert_w_input));
        if u_min < u_0 && u_0 < u_max
            search_points = [search_points; u_0];
        end
    elseif -1/exp(1) <= lambert_w_input && lambert_w_input < 0
        u_0 = real(0.5*value - 0.25*lambertw(lambert_w_input));
        if u_min < u_0 && u_0 < u_max
            search_points = [search_points; u_0];
        end
        u_neg1 = real(0.5*value - 0.25*lambertw(-1, lambert_w_input)); %larger than u_1 because 
            %lambertw(-1, ...) more negative
        if u_min < u_neg1 && u_neg1 < u_max
            search_points = [search_points; u_neg1];
        end
    end
    search_points = [search_points; u_max];
    
    %find the zeros between each search point if there's a sign change
    q_at_search_points = q(search_points, b, d);
    sign_q_at_search_points = sign(q_at_search_points);
    diff_sign_q_at_search_points = diff(sign_q_at_search_points);
    u_candidates = [u_min];
    
    num_tests = 30;
    for search_point_index = 1 : length(search_points) - 1
        count = 0; %number of test points we've tried;
        if diff_sign_q_at_search_points(search_point_index) ~= 0 %if q has different sign from 
            %next search point
            search_start = search_points(search_point_index);
            u_left = search_start;
            q_u_left = q_at_search_points(search_point_index);
            
            search_end = search_points(search_point_index + 1);
            u_right = search_end;
            q_u_right = q_at_search_points(search_point_index);

            u_middle = (u_left + u_right)/2;
            q_u_middle = q(u_middle, b, d);
                
            %Check to see whether the midpoint of u_left and u_right has sign of u_left
            while count < num_tests
                if sign(q_u_middle) == sign(q_u_left) %it's same sign as left, but we want different sign--make midpoint the new left endpoint
                    u_left = u_middle;
                    q_u_left = q_u_middle;
                    count = count + 1;
                else %it's same sign as right, but we want different sign--make midpoint the new right endpoint
                    u_right = u_middle;
                    q_u_right = q_u_middle;
                    count = count + 1;
                end
                u_middle = (u_left + u_right)/2;
                q_u_middle = q(u_middle, b, d);
            end
            u_candidates = [u_candidates; u_middle];
        end
    end
    u_candidates = [u_candidates; u_max];
    
    x_candidates = u_candidates + b;
    function_at_x_candidates = (a*tanh(x_candidates - b) + c)./x_candidates;
    [~, maximizing_x_index] = max(function_at_x_candidates);
    maximizing_x = x_candidates(maximizing_x_index);
    [~, minimizing_x_index] = min(function_at_x_candidates); %"minimizing" is short for "minimizing"
    minimizing_x = x_candidates(minimizing_x_index);
    
    g_abe_max = (a*tanh(maximizing_x - b) + c)./maximizing_x;
    g_abe_min = (a*tanh(minimizing_x - b) + c)./minimizing_x;
end

function output = q(u, b, d)
    output = (1 - d) + (4*u + 4*b - 2*d).*exp(2*u) - (1 + d)*exp(4*u);
end
