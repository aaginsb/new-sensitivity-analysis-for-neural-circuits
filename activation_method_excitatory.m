function y = activation_method_excitatory(x, alpha, beta) 
%ACTIVATION_METHOD 
%Inputs:
    % x = voltage (or Mag) values 
    % max_val = maximum value of the output. 
    % alpha = flatness -- the larger the alpha, the less steep its graph at
        % each point
    % beta = center -- value of input for which output attains half its
        % max. The function's graph is rotationally symmetric about (x,y) = (beta, 0.5*max_val)
%Outputs: y = max firing rate if x is voltage

y = 50 * 0.5 *(1 + tanh((x - beta)/alpha));

end


