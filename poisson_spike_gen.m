function [spike_matrix, t_vec] = poisson_spike_gen(firing_rate, t_sim, n_fibers)
% early drafts of this code was written by Kevin Hannay, and later drafts were written by Alex Ginsberg

%inputs:
    % firing_rate = firing rate (in mHz) of incoming signals
    % t_sim = total time over which we wish to generate a spike train with the specified frequency (ms)
    % n_fibers = number of fibers over which we are simulating input
% outputs:
    % spike_matrix: the ith row corresponds to the ith fiber
        % the jth column corresponds to the jth time bin
        % there is a one in the (i, j) entry of spike_mat if fiber i fires in time bin j
    % t_vec: vector of left endpoints of time bins

% in ms
dt = .1; % time bin (ms)
n_bins = ceil(t_sim/dt);  %number of bins we need for the simulation 
spike_matrix = rand(n_fibers, n_bins) < firing_rate*dt; 
%rand (m, n) outputs a matrix with values sampled uniformly from [0, 1]
    % rand(m, n) < x outputs an m-by-n matrix with values sampled uniformly from
    % [0, 1] and then replaced with a 0 if the value is larger than x and 1
    % if the value is smaller than x. A 1 in the (i, j) entry of
    % spike_matrix represents a spike in the ith fiber during the jth time
    % interval. If 0 < firing_rate*dt < 1, there is a 1 in any particular
    % entry with probability firing_rate*dt, and the expected number of
    % spikes in millisecond is then n_bins*firing_rate*dt =
    % 1/dt*firing_rate*dt = firing_rate (mHZ). Letting dt tend to 0 yields
    % a poisson process in each row

t_vec = 0: dt : t_sim - dt; % vector of left endpoints of time bins

end